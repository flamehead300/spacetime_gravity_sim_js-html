#!/usr/bin/env python3
"""
Gravity Simulator – Python-authoritative physics, browser renders.
Zero pip dependencies.  Run with:  python gravity_web.py
"""
# Import JSON for web serialization, math for physics equations, random for body generation
import json, math, random, threading, webbrowser
# Import basic HTTP server components to serve the physics state to a web frontend
from http.server import HTTPServer, BaseHTTPRequestHandler

# Define the port the local web server will run on
PORT = 8765

# ═══════════════════════════════════════════════════════════════════
# Physics engine  (ported 1-to-1 from the browser JS)
# ═══════════════════════════════════════════════════════════════════
# Define the minimum and maximum allowable mass for a generated body
MASS_MIN, MASS_MAX = 1, 1_000_000
# Define the minimum and maximum visual/physical radius of a body
R_MIN, R_MAX = 1.0, 10.0
# Pre-calculate cube roots of mass limits for volume-to-radius scaling
_CBRT_MIN = MASS_MIN ** (1/3)
_CBRT_MAX = MASS_MAX ** (1/3)
# Reference delta-time used to normalize velocity damping (friction)
_DAMPING_REF_DT = 1 / 360
# Maximum number of collision resolution iterations per physics step to prevent infinite loops
_MAX_COLL = 8

def _default_mass_to_r(m):
    # Clamp the input mass 'm' between the defined minimum and maximum mass
    c = max(MASS_MIN, min(MASS_MAX, m))
    # Interpolate the radius based on the cube root of the mass (since volume ∝ mass ∝ r^3)
    return R_MIN + (c ** (1/3) - _CBRT_MIN) * (R_MAX - R_MIN) / (_CBRT_MAX - _CBRT_MIN)

def _rs(G, m, c):
    # Calculate the Schwarzschild radius (event horizon) formula: 2Gm / c^2
    return (2 * G * m) / (c * c) if c > 0 else 0

def _rand_color():
    # Generate a random RGB color string, keeping values above 140 to ensure bright, visible colors
    return "rgb({},{},{})".format(
        random.randint(140, 255), random.randint(140, 255), random.randint(140, 255))


class Sim:
    def __init__(self):
        # List to hold the state dictionaries of all active physics bodies
        self.bodies: list[dict] = []
        # Monotonically increasing counter to assign unique IDs to each body
        self.next_id = 0
        # Gravitational constant (scaled up for visual simulation purposes)
        self.G            = 10.0
        # Softening parameter to prevent infinite forces when distance approaches zero
        self.softening    = 1.5
        # Physics time-step per frame (dt)
        self.dt           = 0.002
        # Velocity multiplier applied per step to simulate vacuum drag/friction (1.0 = no drag)
        self.damping      = 1.0
        # Absolute hard cap on velocity to prevent physics explosions/tunneling
        self.max_speed    = 8000.0
        # Speed of light in the simulation, used for relativistic limits and 1PN gravity
        self.c            = 8000.0
        # Boolean flag to toggle simple relativistic velocity capping
        self.use_rel      = True
        # Boolean flag to toggle First Post-Newtonian (1PN) general relativity gravity simulation
        self.use_1pn      = False
        # Optional higher-order correction layered on top of the existing 1PN approximation
        self.use_2pn      = False
        # Enable simple gravitational-wave radiation reaction for sufficiently massive pairs
        self.gw_damping   = True
        # Only apply GW losses to heavy-body proxies so low-mass particle clouds stay responsive
        self.gw_mass_threshold = 100_000.0
        # Tidal breakup of the lighter body when it enters the Roche limit of a heavier one
        self.tidal_frag   = True
        # Adapt the step size downward when the system gets dynamically stiff
        self.adaptive_dt  = True
        self.dt_min       = 0.0005
        self.dt_max       = 0.01
        # Optional alternate display scaling for rocky / gas / degenerate bodies
        self.realistic_radii = False
        # Optional symplectic-style integrator; disabled by default because it doubles force evaluations
        self.use_leapfrog = False
        # Spawn a small ejecta fragment from energetic inelastic mergers
        self.ejecta_enabled = False
        # Optional background halo potential
        self.dark_matter  = False
        self.dm_scale     = 1000.0
        self.dm_mass      = 1e7
        # Coefficient of restitution (0.0 = perfectly inelastic/merging, 1.0 = perfectly elastic/bouncing)
        self.restitution  = 0.0
        # Boolean flag to allow bodies to break apart upon high-speed collisions
        self.frag_enabled = False
        # Kinetic energy/velocity threshold required to trigger a fragmentation event
        self.frag_threshold = 250.0
        # Skip collision resolution once the scene is large enough that all-pairs overlap tests dominate.
        self.collision_limit = 250
        # Hard cap keeps the demo responsive instead of allowing the scene to spiral into unusable counts.
        self.max_bodies = 1200
        # Threading lock to prevent the HTTP server from reading state while the physics thread is mutating it
        self.lock = threading.Lock()

    def _mass_to_r(self, m):
        c = max(MASS_MIN, min(MASS_MAX, float(m)))
        if not self.realistic_radii:
            return _default_mass_to_r(c)

        rocky_max = 100_000.0
        gas_max = 1_000_000.0
        if c <= rocky_max:
            return max(R_MIN, min(R_MAX, c ** (1/3)))
        if c <= gas_max:
            rocky_r = rocky_max ** (1/3)
            gas_r = rocky_r * (c / rocky_max) ** 0.55
            return max(R_MIN, min(R_MAX, gas_r))

        degen_r = 10.0 * (gas_max / c) ** (1/3)
        return max(R_MIN, min(R_MAX, degen_r))

    def _refresh_body_radii(self):
        for b in self.bodies:
            b['r'] = self._mass_to_r(b['m'])

    def _estimate_dt(self, bs):
        dt = max(self.dt_min, min(self.dt_max, self.dt))
        if not self.adaptive_dt or len(bs) < 2:
            return dt

        soft2 = max(self.softening * self.softening, 1e-9)
        total_mass = sum(b['m'] for b in bs)
        max_a = 0.0
        for b in bs:
            a_est = self.G * max(total_mass - b['m'], 0.0) / soft2
            if a_est > max_a:
                max_a = a_est
        if max_a <= 0:
            return dt

        dt_est = 0.2 * math.sqrt(soft2 / max_a)
        return max(self.dt_min, min(self.dt_max, min(dt, dt_est)))

    def _dm_acceleration(self, x, y):
        r2 = x * x + y * y
        if r2 < 1e-9:
            return 0.0, 0.0
        scale2 = self.dm_scale * self.dm_scale
        # Plummer-like halo keeps the central force bounded while still acting as a background potential.
        factor = -self.G * self.dm_mass / ((r2 + scale2) ** 1.5)
        return x * factor, y * factor

    def _tidal_fragment(self, i, j):
        bi, bj = self.bodies[i], self.bodies[j]
        if bi['m'] > bj['m']:
            i, j = j, i
            bi, bj = self.bodies[i], self.bodies[j]

        if bi['m'] <= 0:
            return False

        r_roche = bi['r'] * (2 * bj['m'] / bi['m']) ** (1.0 / 3.0)
        dx = bj['x'] - bi['x']
        dy = bj['y'] - bi['y']
        dist = math.hypot(dx, dy)
        if dist <= 0 or dist >= r_roche:
            return False

        tidally_disrupted = self.bodies.pop(i)
        frag_m = max(1.0, tidally_disrupted['m'] / 3.0)
        fr = self._mass_to_r(frag_m)
        tangent = math.atan2(dy, dx) + math.pi / 2
        rel_speed = math.hypot(
            tidally_disrupted['vx'] - bj['vx'],
            tidally_disrupted['vy'] - bj['vy'],
        )
        spread = max(rel_speed, self.frag_threshold * 0.15)
        for k in range(3):
            angle = tangent + (k - 1) * 0.65 + (random.random() - 0.5) * 0.2
            speed = spread * (0.6 + random.random() * 0.4)
            self.add_body(
                tidally_disrupted['x'] + math.cos(angle) * (fr + 1) * 2.0,
                tidally_disrupted['y'] + math.sin(angle) * (fr + 1) * 2.0,
                tidally_disrupted['vx'] + math.cos(angle) * speed,
                tidally_disrupted['vy'] + math.sin(angle) * speed,
                frag_m,
            )
        return True

    # ── body management ─────────────────────────────────────────────
    def add_body(self, x, y, vx, vy, m):
        if len(self.bodies) >= self.max_bodies:
            return False
        # Create a dictionary representing the physical and visual properties of the new body
        b = dict(id=self.next_id, x=x, y=y, vx=vx, vy=vy,
                 m=float(m), r=self._mass_to_r(m), color=_rand_color())
        # Increment the unique ID counter
        self.next_id += 1
        # Append the new body to the active simulation array
        self.bodies.append(b)
        return True

    def _merge(self, si, di):
        # Retrieve the source (s) and destination (d) bodies from the array
        s, d = self.bodies[si], self.bodies[di]
        # Extract their respective masses
        m1, m2 = s['m'], d['m']
        # Calculate the total mass (ensure it's not zero to avoid division by zero)
        mt = m1 + m2 or 1
        vrelx = s['vx'] - d['vx']
        vrely = s['vy'] - d['vy']
        rel_ke = 0.5 * (m1 * m2 / mt) * (vrelx * vrelx + vrely * vrely)
        ejecta_mass = 0.0
        if self.ejecta_enabled and self.restitution < 0.05 and mt > 2.0 and rel_ke > self.frag_threshold * 0.5:
            ejecta_mass = min(max(1.0, mt * 0.01), mt * 0.25)
        merged_mass = mt - ejecta_mass
        # Calculate the new X position using the center of mass formula
        d['x']  = (m1*s['x']  + m2*d['x'])  / mt
        # Calculate the new Y position using the center of mass formula
        d['y']  = (m1*s['y']  + m2*d['y'])  / mt
        # Calculate the new X velocity (Conservation of Momentum: m1v1 + m2v2 = mtvt)
        d['vx'] = (m1*s['vx'] + m2*d['vx']) / mt
        # Calculate the new Y velocity
        d['vy'] = (m1*s['vy'] + m2*d['vy']) / mt
        # Set the destination body's mass to the combined total mass
        d['m']  = merged_mass
        # Recalculate the destination body's physical radius based on its new mass
        d['r']  = self._mass_to_r(merged_mass)
        # Remove the source body from the simulation since it has been absorbed
        self.bodies.pop(si)
        if ejecta_mass > 0:
            angle = random.uniform(0, 2 * math.pi)
            speed = math.sqrt(2 * rel_ke / ejecta_mass) * 0.25
            self.add_body(
                d['x'],
                d['y'],
                d['vx'] + math.cos(angle) * speed,
                d['vy'] + math.sin(angle) * speed,
                ejecta_mass,
            )

    def _fragment(self, si, di):
        # Retrieve the two colliding bodies
        s, d = self.bodies[si], self.bodies[di]
        # Extract their masses
        m1, m2 = s['m'], d['m']
        # Calculate total mass of the system
        mt = m1 + m2
        # Calculate the X coordinate of the center of mass
        cmx  = (m1*s['x']  + m2*d['x'])  / mt
        # Calculate the Y coordinate of the center of mass
        cmy  = (m1*s['y']  + m2*d['y'])  / mt
        # Calculate the X velocity of the center of mass
        cmvx = (m1*s['vx'] + m2*d['vx']) / mt
        # Calculate the Y velocity of the center of mass
        cmvy = (m1*s['vy'] + m2*d['vy']) / mt
        # Calculate relative velocity components between the two bodies
        dvx  = s['vx'] - d['vx'];  dvy = s['vy'] - d['vy']
        # Calculate the relative kinetic energy using the reduced mass (m1*m2/mt)
        rel_ke   = 0.5 * (m1*m2/mt) * (dvx*dvx + dvy*dvy)
        # Determine ejection speed of fragments based on relative KE and a baseline threshold
        ej_speed = math.sqrt(max(0, 2*rel_ke/mt)) * 0.5 + self.frag_threshold * 0.15
        # Identify the higher and lower indices to safely pop elements from the list without shifting issues
        hi, lo = max(si, di), min(si, di)
        # Remove both original colliding bodies from the simulation
        self.bodies.pop(hi);  self.bodies.pop(lo)
        # Distribute the total mass equally into 3 fragments (minimum mass 1.0)
        frag_m = max(1.0, mt / 3)
        # Calculate the radius of these new fragments
        fr     = self._mass_to_r(frag_m)
        # Loop 3 times to create the 3 fragment bodies
        for k in range(3):
            # Calculate an ejection angle separated by roughly 120 degrees (2pi/3), with slight randomization
            angle = 2*math.pi*k/3 + (random.random()-0.5)*0.4
            # Randomize the ejection speed slightly for varied scatter patterns
            spd   = ej_speed * (0.8 + random.random()*0.4)
            # Add the new fragment body to the simulation, offsetting position so they don't instantly collide
            self.add_body(
                cmx + math.cos(angle)*(fr+1)*2.5, # Offset X
                cmy + math.sin(angle)*(fr+1)*2.5, # Offset Y
                cmvx + math.cos(angle)*spd,       # Add ejection velocity X to center of mass velocity
                cmvy + math.sin(angle)*spd,       # Add ejection velocity Y to center of mass velocity
                frag_m)                           # Assign fragment mass

    # ── collision resolution ─────────────────────────────────────────
    def _collect_pairs(self):
        # Cache the bodies list locally for faster access
        bs = self.bodies
        # Get the total number of bodies
        n  = len(bs)
        # Initialize an empty list to store colliding pairs
        out = []
        # Outer loop iterates through every body (except the last)
        for i in range(n-1):
            bi = bs[i]
            # Inner loop iterates through all subsequent bodies to check for unique pairs
            for j in range(i+1, n):
                bj = bs[j]
                # Calculate delta X and delta Y between the two bodies
                dx = bj['x']-bi['x'];  dy = bj['y']-bi['y']
                # Calculate the minimum distance required to NOT be colliding (sum of radii)
                lim = bi['r']+bj['r'];  
                # Calculate squared distance to avoid expensive square root unless necessary
                d2 = dx*dx+dy*dy
                # Broad-phase check: If squared distance >= squared limit, they aren't colliding. Skip.
                if d2 >= lim*lim: continue
                # Narrow-phase: Calculate actual Euclidean distance
                dist = math.sqrt(d2)
                # Append overlap amount, unique pair hash, indices, distance, and deltas
                out.append((lim-dist, bi['id']*2_000_000+bj['id'], i, j, dist, dx, dy))
        # Sort collisions descending by overlap amount to resolve the worst collisions first
        out.sort(reverse=True)   # id-pair tie-break determinism ensures identical cross-platform behavior
        return out

    def _resolve(self):
        # Iterate up to maximum allowed collision resolutions per frame
        for _ in range(_MAX_COLL):
            # Get the current list of colliding pairs
            pairs = self._collect_pairs()
            # If no bodies are overlapping, exit the resolution loop early
            if not pairs: break
            # Extract the data for the most overlapping pair (index 0)
            _, _, i, j, dist, dx, dy = pairs[0]
            # Retrieve the actual body dictionaries
            bi, bj = self.bodies[i], self.bodies[j]
            # Ensure distance is not zero to prevent division by zero in normal calculation
            d  = max(dist, 1e-9)
            # Calculate collision normal vector components (normalized direction of collision)
            nx, ny = dx/d, dy/d
            # Calculate relative velocity components
            dvx = bi['vx']-bj['vx'];  dvy = bi['vy']-bj['vy']
            # Calculate relative velocity along the normal (dot product)
            dvn = dvx*nx + dvy*ny

            if self.tidal_frag and self.restitution < 0.05:
                if self._tidal_fragment(i, j):
                    continue

            # If restitution is very low, treat it as an inelastic collision (merge or fragment)
            if self.restitution < 0.05:
                # If fragmentation is enabled and impact speed exceeds the threshold, break them apart
                if self.frag_enabled and abs(dvn) > self.frag_threshold:
                    self._fragment(i, j)
                # Otherwise, merge the smaller body into the larger body
                else:
                    # Determine which index has the larger mass to act as the destination
                    dst = i if bi['m'] >= bj['m'] else j
                    # The other index acts as the source being absorbed
                    src = j if dst == i else i
                    # Perform the merge operation
                    self._merge(src, dst)
            # If restitution is high, treat it as an elastic collision (bounce)
            else:
                # Only apply bounce impulse if bodies are moving towards each other (dvn > 0)
                if dvn > 0:
                    mi, mj = bi['m'], bj['m']
                    # Calculate impulse scalar (j) based on restitution, relative velocity, and masses
                    ji = (1+self.restitution)*dvn*mi*mj/(mi+mj)
                    # Apply impulse to body 'i', altering its velocity away from 'j'
                    bi['vx'] -= ji/mi*nx;  bi['vy'] -= ji/mi*ny
                    # Apply equal and opposite impulse to body 'j'
                    bj['vx'] += ji/mj*nx;  bj['vy'] += ji/mj*ny
                # Recalculate required separation distance (sum of radii)
                lim_now = bi['r']+bj['r']
                # Calculate how deeply they intersect, add a small buffer (0.05) to push them fully apart
                sep = lim_now - dist + 0.05
                # Total mass for positional correction weighting
                mt  = bi['m']+bj['m']
                # Push body 'i' out of the collision based on its mass ratio (heavier moves less)
                bi['x'] -= nx*sep*bj['m']/mt;  bi['y'] -= ny*sep*bj['m']/mt
                # Push body 'j' out of the collision in the opposite direction
                bj['x'] += nx*sep*bi['m']/mt;  bj['y'] += ny*sep*bi['m']/mt

    # ── integrator ───────────────────────────────────────────────────
    def _compute_accelerations(self, bs, dt):
        n = len(bs)
        G = self.G
        c = self.c
        eps2 = self.softening ** 2
        ax = [0.0] * n
        ay = [0.0] * n

        for i in range(n - 1):
            bi = bs[i]
            for j in range(i + 1, n):
                bj = bs[j]
                dx = bj['x'] - bi['x']
                dy = bj['y'] - bi['y']
                rs2 = dx * dx + dy * dy + eps2
                if rs2 < 1e-9:
                    continue

                if not self.use_rel:
                    ir3 = rs2 ** -1.5
                    si = G * bj['m'] * ir3
                    sj = G * bi['m'] * ir3
                    ax[i] += dx * si
                    ay[i] += dy * si
                    ax[j] -= dx * sj
                    ay[j] -= dy * sj

                elif not self.use_1pn:
                    rs = math.sqrt(rs2)
                    di = max(rs - _rs(G, bj['m'], c), 1)
                    dj = max(rs - _rs(G, bi['m'], c), 1)
                    ir = 1 / rs
                    si = G * bj['m'] * ir / (di * di)
                    sj = G * bi['m'] * ir / (dj * dj)
                    ax[i] += dx * si
                    ay[i] += dy * si
                    ax[j] -= dx * sj
                    ay[j] -= dy * sj

                else:
                    rs = math.sqrt(rs2)
                    ir = 1.0 / rs
                    ir3 = ir / rs2
                    si = G * bj['m'] * ir3
                    sj = G * bi['m'] * ir3
                    ax[i] += dx * si
                    ay[i] += dy * si
                    ax[j] -= dx * sj
                    ay[j] -= dy * sj

                    c2 = c * c
                    nx = dx * ir
                    ny = dy * ir
                    vi2 = bi['vx'] ** 2 + bi['vy'] ** 2
                    vj2 = bj['vx'] ** 2 + bj['vy'] ** 2
                    vivj = bi['vx'] * bj['vx'] + bi['vy'] * bj['vy']
                    vi_n = bi['vx'] * nx + bi['vy'] * ny
                    vj_n = bj['vx'] * nx + bj['vy'] * ny
                    Gmi_r = G * bi['m'] * ir
                    Gmj_r = G * bj['m'] * ir
                    ir2 = ir * ir
                    dvx2 = bi['vx'] - bj['vx']
                    dvy2 = bi['vy'] - bj['vy']

                    fi_f = G * bj['m'] * ir2 / c2
                    fi_s = -vi2 - 2 * vj2 + 4 * vivj + 1.5 * vj_n ** 2 + 5 * Gmi_r + 4 * Gmj_r
                    fi_v = 4 * vi_n - 3 * vj_n
                    pn1_ix = fi_f * (nx * fi_s + dvx2 * fi_v)
                    pn1_iy = fi_f * (ny * fi_s + dvy2 * fi_v)
                    ax[i] += pn1_ix
                    ay[i] += pn1_iy

                    fj_f = G * bi['m'] * ir2 / c2
                    fj_s = -vj2 - 2 * vi2 + 4 * vivj + 1.5 * vi_n ** 2 + 5 * Gmj_r + 4 * Gmi_r
                    fj_v = 4 * vj_n - 3 * vi_n
                    pn1_jx = fj_f * (-nx * fj_s + dvx2 * fj_v)
                    pn1_jy = fj_f * (-ny * fj_s + dvy2 * fj_v)
                    ax[j] += pn1_jx
                    ay[j] += pn1_jy

                    if self.use_2pn:
                        # This remains an educational approximation rather than a full 2PN EIH expansion:
                        # scale the already-computed 1PN correction by an extra compactness/speed factor.
                        pn_boost = min(0.35, 2.0 * (0.5 * (vi2 + vj2) + Gmi_r + Gmj_r) / c2)
                        ax[i] += pn_boost * pn1_ix
                        ay[i] += pn_boost * pn1_iy
                        ax[j] += pn_boost * pn1_jx
                        ay[j] += pn_boost * pn1_jy

                if self.use_rel and self.gw_damping and (
                    bi['m'] >= self.gw_mass_threshold or bj['m'] >= self.gw_mass_threshold
                ):
                    r = math.sqrt(rs2)
                    if r > 0:
                        vrelx = bj['vx'] - bi['vx']
                        vrely = bj['vy'] - bi['vy']
                        vrel2 = vrelx * vrelx + vrely * vrely
                        if vrel2 > 1e-12:
                            vrel = math.sqrt(vrel2)
                            ux = vrelx / vrel
                            uy = vrely / vrel
                            beta2 = min(1.0, vrel2 / max(c * c, 1e-9))
                            mi = max(bi['m'], 1e-9)
                            mj = max(bj['m'], 1e-9)
                            dE_dt = (32.0 / 5.0) * (
                                G ** 4 * mi * mi * mj * mj * (mi + mj)
                                / max(c ** 5 * r ** 5, 1e-9)
                            )
                            force = dE_dt / vrel
                            force *= beta2
                            reduced_mass = (mi * mj) / (mi + mj)
                            max_force = 0.25 * reduced_mass * vrel / max(dt, 1e-9)
                            force = min(force, max_force)
                            ax[i] += ux * (force / mi)
                            ay[i] += uy * (force / mi)
                            ax[j] -= ux * (force / mj)
                            ay[j] -= uy * (force / mj)

        if self.dark_matter:
            for i, b in enumerate(bs):
                dax, day = self._dm_acceleration(b['x'], b['y'])
                ax[i] += dax
                ay[i] += day

        return ax, ay

    def step(self, n=1):
        # Step the physics engine forward 'n' times
        for _ in range(n):
            self._step_once()

    def _step_once(self):
        bs = self.bodies
        if not bs:
            return

        do_collisions = len(bs) <= self.collision_limit
        if do_collisions:
            self._resolve()
        dt = self._estimate_dt(bs)
        drag = self.damping ** (dt / _DAMPING_REF_DT)
        vcap = min(self.max_speed, self.c) if self.use_rel else self.max_speed
        vcap2 = vcap * vcap

        if self.use_leapfrog:
            ax, ay = self._compute_accelerations(bs, dt)
            for i, b in enumerate(bs):
                b['vx'] += ax[i] * dt * 0.5
                b['vy'] += ay[i] * dt * 0.5
            for b in bs:
                b['x'] += b['vx'] * dt
                b['y'] += b['vy'] * dt
            ax, ay = self._compute_accelerations(bs, dt)
            for i, b in enumerate(bs):
                b['vx'] = (b['vx'] + ax[i] * dt * 0.5) * drag
                b['vy'] = (b['vy'] + ay[i] * dt * 0.5) * drag
                v2 = b['vx'] * b['vx'] + b['vy'] * b['vy']
                if v2 > vcap2:
                    s = vcap / math.sqrt(v2)
                    b['vx'] *= s
                    b['vy'] *= s
        else:
            ax, ay = self._compute_accelerations(bs, dt)
            for i, b in enumerate(bs):
                b['vx'] = (b['vx'] + ax[i] * dt) * drag
                b['vy'] = (b['vy'] + ay[i] * dt) * drag
                v2 = b['vx'] * b['vx'] + b['vy'] * b['vy']
                if v2 > vcap2:
                    s = vcap / math.sqrt(v2)
                    b['vx'] *= s
                    b['vy'] *= s
                b['x'] += b['vx'] * dt
                b['y'] += b['vy'] * dt

        if do_collisions:
            self._resolve()
        return


    # ── serialisation ────────────────────────────────────────────────
    def snapshot(self):
        # Return a complete dictionary representation of the simulation state (ideal for JSON.dumps to a web client)
        return {
            # Map over bodies to create a clean list of dictionaries
            'bodies': [
                {'id':b['id'],'x':b['x'],'y':b['y'],'vx':b['vx'],'vy':b['vy'],
                 'm':b['m'],'r':b['r'],'color':b['color']}
                for b in self.bodies],
            # Pass out current engine parameters so the frontend UI stays synchronized
            'params': {
                'G':self.G,'softening':self.softening,'dt':self.dt,
                'damping':self.damping,'max_speed':self.max_speed,'c':self.c,
                'use_rel':self.use_rel,'use_1pn':self.use_1pn,'use_2pn':self.use_2pn,
                'gw_damping':self.gw_damping,'tidal_frag':self.tidal_frag,
                'adaptive_dt':self.adaptive_dt,'dt_min':self.dt_min,'dt_max':self.dt_max,
                'realistic_radii':self.realistic_radii,'use_leapfrog':self.use_leapfrog,
                'ejecta_enabled':self.ejecta_enabled,
                'dark_matter':self.dark_matter,'dm_scale':self.dm_scale,'dm_mass':self.dm_mass,
                'restitution':self.restitution,
                'frag_enabled':self.frag_enabled,'frag_threshold':self.frag_threshold,
            }
        }

    def set_params(self, d: dict):
        # Define a set of keys that should safely be cast to floats
        float_keys = {
            'G','softening','dt','damping','max_speed','c',
            'frag_threshold','restitution','dt_min','dt_max','dm_scale','dm_mass'
        }
        # Define a set of keys that should safely be cast to booleans
        bool_keys  = {
            'use_rel','use_1pn','use_2pn','gw_damping','tidal_frag',
            'adaptive_dt','realistic_radii','use_leapfrog','ejecta_enabled',
            'dark_matter','frag_enabled'
        }
        # Iterate over incoming dictionary payload (typically from an HTTP POST request)
        refresh_radii = False
        for k, v in d.items():
            # If key matches float list, dynamically set the attribute on the Sim instance
            if k in float_keys: setattr(self, k, float(v))
            # If key matches boolean list, dynamically set the attribute
            elif k in bool_keys: setattr(self, k, bool(v))
            if k == 'realistic_radii':
                refresh_radii = True
        if self.dt_min > self.dt_max:
            self.dt_min, self.dt_max = self.dt_max, self.dt_min
        if self.use_2pn:
            self.use_1pn = True
        elif not self.use_1pn:
            self.use_2pn = False
        if refresh_radii:
            self._refresh_body_radii()

    def clear(self):
        # Instantly wipe all bodies from the simulation array
        self.bodies.clear()
# ═══════════════════════════════════════════════════════════════════ # Top visual divider for the HTTP handler section
# HTTP handler                                                        # Section title indicating backend server logic
# ═══════════════════════════════════════════════════════════════════ # Bottom visual divider for the HTTP handler section
SIM = Sim() # Instantiate the global physics simulation engine object

class _H(BaseHTTPRequestHandler): # Define a custom request handler inheriting from Python's built-in HTTP server
    def _read_json(self): # Helper method to extract and parse JSON from incoming POST requests
        n = int(self.headers.get('Content-Length', 0)) # Retrieve the Content-Length header, defaulting to 0 if missing
        return json.loads(self.rfile.read(n)) if n else {} # Read 'n' bytes from the request body and parse JSON, or return empty dict

    def _send(self, data: dict, status=200): # Helper method to send a JSON payload back to the client
        body = json.dumps(data).encode() # Serialize the Python dictionary into a JSON string, then encode it into bytes
        self.send_response(status) # Send the HTTP status code (defaults to 200 OK)
        self.send_header('Content-Type', 'application/json') # Declare the response content type as JSON
        self.send_header('Content-Length', str(len(body))) # Declare the exact byte length of the response payload
        self.end_headers() # Send the blank line indicating the end of HTTP headers
        self.wfile.write(body) # Write the encoded JSON bytes directly to the output stream

    def do_GET(self): # Override method to handle all incoming HTTP GET requests
        if self.path == '/state': # Check if the client is requesting the current simulation state
            with SIM.lock: # Acquire the thread lock to safely read data from the concurrently running simulation thread
                self._send(SIM.snapshot()) # Serialize and send a snapshot of the current simulation variables
        else: # Fallback for all other GET routes (primarily the root '/' to load the UI)
            html = HTML.encode() # Encode the global raw HTML string into bytes
            self.send_response(200) # Send a 200 OK HTTP status code
            self.send_header('Content-Type', 'text/html; charset=utf-8') # Declare the response as UTF-8 encoded HTML
            self.send_header('Content-Length', str(len(html))) # Provide the exact byte length of the HTML string
            self.end_headers() # End the HTTP headers section
            self.wfile.write(html) # Write the HTML bytes to the stream, rendering the webpage in the browser

    def do_POST(self): # Override method to handle all incoming HTTP POST requests (used for client-to-server commands)
        p = self.path.split('?')[0] # Extract the base endpoint path, discarding any URL query parameters
        d = self._read_json() # Call the helper to parse the JSON request body into a dictionary
        with SIM.lock: # Acquire the thread lock to prevent data corruption while mutating simulation state
            if p == '/step': # Endpoint: Advance the simulation forward manually
                SIM.step(int(d.get('n', 1))) # Step the simulation by 'n' ticks (defaulting to 1 if not provided)
                self._send(SIM.snapshot()) # Return the updated simulation state to the client
            elif p == '/spawn': # Endpoint: Spawn a single celestial body
                SIM.add_body(d['x'], d['y'], d['vx'], d['vy'], d['m']) # Add a body using the provided x/y coords, velocities, and mass
                self._send(SIM.snapshot()) # Return the updated state confirming the body was added
            elif p == '/spawn_batch': # Endpoint: Spawn multiple bodies in a single request (optimization)
                added = 0 # Initialize a counter to track successfully added bodies
                for b in d.get('bodies', []): # Iterate through the 'bodies' array in the JSON payload
                    if SIM.add_body(b['x'], b['y'], b['vx'], b['vy'], b['m']): # Add each individual body to the simulation
                        added += 1 # Increment the tracking counter
                self._send({'ok': True, 'added': added, 'count': len(SIM.bodies)}) # Return success confirmation and total body count
            elif p == '/clear': # Endpoint: Erase all bodies from the simulation
                SIM.clear() # Execute the clear command on the engine
                self._send(SIM.snapshot()) # Return the new (empty) state
            elif p == '/set_params': # Endpoint: Modify global physics parameters
                SIM.set_params(d) # Pass the dictionary of new parameters to the simulation engine
                self._send({'ok': True, 'params': SIM.snapshot()['params']}) # Confirm success and return the active parameters
            else: # Fallback: The requested POST endpoint does not exist
                self._send({'error': 'unknown route'}, 404) # Send a 404 Not Found error wrapped in JSON

    def log_message(self, *_): pass # Silence the default BaseHTTPRequestHandler logging to keep the terminal output clean


# ═══════════════════════════════════════════════════════════════════ # Top visual divider for the Client UI section
# Browser client                                                      # Section title indicating Frontend code
# ═══════════════════════════════════════════════════════════════════ # Bottom visual divider for the Client UI section
HTML = r"""<!doctype html> <!-- Declare standard HTML5 document type -->
<html lang="en"> <!-- Start the HTML document, setting language to English -->
<head> <!-- Open the document head for metadata and styling -->
<meta charset="utf-8"> <!-- Specify UTF-8 character encoding -->
<meta name="viewport" content="width=device-width,initial-scale=1"> <!-- Ensure responsive scaling on mobile devices -->
<title>STELLAR DYNAMICS // GRAVITY ENGINE</title> <!-- Set the text displayed in the browser tab -->
<style> /* Open inline CSS style block */
*{box-sizing:border-box;margin:0;padding:0} /* Universal reset: Use border-box sizing and remove default margins/padding */
html,body{height:100%;overflow:hidden;background:#050a10} /* Make body fill viewport height, hide scrollbars, set dark space-blue background */
body{display:flex;font-family:'Courier New',Courier,monospace;color:#00e5cc} /* Flex layout, monospace font for UI, cyan text color */
#cw{flex:1;min-width:0;padding:8px 0 8px 8px} /* Canvas wrapper: Grow to fill remaining space, add padding */
canvas{display:block;width:100%;height:100%;outline:1px solid #003d33;cursor:crosshair} /* Canvas: Fill wrapper, add dark green border, use crosshair cursor */
#panel{ /* Sidebar control panel styles */
  width:228px;flex-shrink:0;background:#070e16;border-left:1px solid #00e5cc; /* Fixed width of 228px, don't shrink, dark background, cyan left border */
  display:flex;flex-direction:column;overflow-y:auto;padding:10px 8px 14px; /* Stack children vertically, enable vertical scroll, set padding */
  scrollbar-width:thin;scrollbar-color:#003d33 #070e16} /* Firefox scrollbar styling (dark green thumb, dark blue track) */
.hdr{font-size:13px;font-weight:bold;margin-bottom:1px} /* Main header text: 13px, bold, minimal bottom margin */
.sub{font-size:10px;color:#005544;margin-bottom:10px} /* Subheader text: 10px, dark green color, spacing below */
.sec{border:1px solid #00e5cc;padding:7px 7px 3px;margin-bottom:7px;position:relative} /* UI Section box: Cyan border, inner padding, positioned relatively for title overlap */
.sec-t{position:absolute;top:-8px;left:6px;background:#070e16;padding:0 3px;font-size:10px;font-weight:bold} /* Section title: Break border outline, sit on top edge, match background */
.lbl{display:flex;justify-content:space-between;font-size:10px;margin-bottom:1px} /* Label container: Flex to push name left and value right */
.lbl span{color:#00b09c} /* Specific style for value readout spans: slightly darker cyan */
.info{font-size:10px;color:#00b09c;line-height:1.6;margin-bottom:3px} /* Helper text style: small, readable line height, bottom margin */
label.chk{display:flex;align-items:center;gap:5px;font-size:10px;cursor:pointer;margin-bottom:4px} /* Checkbox container: Flex row, center items, add gap between box and text */
input[type=checkbox]{accent-color:#00e5cc;width:12px;height:12px;cursor:pointer} /* Checkbox styling: Cyan checkmark color, fixed square size */
input[type=range]{ /* Universal slider styling */
  -webkit-appearance:none;appearance:none;width:100%;height:4px; /* Strip native OS styling, full width, 4px thick */
  background:#003d33;border-radius:2px;outline:none;margin-bottom:5px;cursor:pointer} /* Dark green track, rounded edges, remove focus outline */
input[type=range]::-webkit-slider-thumb{ /* Webkit (Chrome/Safari) slider thumb styling */
  -webkit-appearance:none;width:10px;height:16px;background:#00e5cc;border-radius:2px;cursor:pointer} /* Strip native, 10x16px cyan rectangle */
input[type=range]::-moz-range-thumb{ /* Firefox slider thumb styling */
  width:10px;height:16px;background:#00e5cc;border-radius:2px;border:none;cursor:pointer} /* Replicate webkit look for cross-browser consistency */
button{ /* Universal button styling */
  display:block;width:100%;background:#070e16;color:#00e5cc;border:1px solid #00e5cc; /* Full width block, dark bg, cyan text and border */
  font-family:inherit;font-size:11px;font-weight:bold;padding:5px 6px; /* Inherit monospace font, 11px bold text, padding */
  margin-bottom:4px;text-align:left;cursor:crosshair} /* Spacing below, left-aligned text, thematic crosshair cursor */
button:hover{background:#003d33} /* Highlight button background to dark green on hover */
button.warn{color:#ffaa00;border-color:#ffaa00} /* 'Warning' button class (used for clear): Orange text and border */
button.warn:hover{background:#3d2600} /* Warning button hover: Dark orange/brown background */
button.active{color:#ffaa00;border-color:#ffaa00} /* 'Active' toggle state: Highlight button in orange */
.tele{font-size:11px;line-height:1.8} /* Telemetry stats text: 11px font, tall line height for readability */
</style> <!-- Close inline CSS block -->
</head> <!-- Close HTML head section -->
<body> <!-- Open HTML body section -->
<div id="cw"><canvas id="sim"></canvas></div> <!-- Render area wrapper and the WebGL/2D Canvas element itself -->
<div id="panel"> <!-- Sidebar controls container -->
  <div class="hdr">&#9656; STELLAR DYNAMICS</div> <!-- App title using Unicode triangle right arrow -->
  <div class="sub">&nbsp;&nbsp;GRAVITY ENGINE v3.1 [JS]</div> <!-- App subtitle indicating browser-side physics -->

  <div class="sec"> <!-- Settings group: Body spawning parameters -->
    <span class="sec-t">BODY PARAMS</span> <!-- Label breaking the top border of the box -->
    <div class="lbl">SPAWN MASS <span id="mass-val">10,000</span></div> <!-- Label indicating mass slider purpose, with dynamic value display -->
    <input type="range" id="mass" min="1" max="1000000" step="1" value="10000"> <!-- Slider for spawning mass: 1 to 1,000,000 -->
    <div class="info" id="body-info"></div> <!-- Empty div reserved for injecting dynamic body info via JS later -->
  </div> <!-- Close Body Params section -->

  <div class="sec"> <!-- Settings group: Physics Engine constants -->
    <span class="sec-t">PHYSICS ENGINE</span> <!-- Label breaking top border -->
    <div class="lbl">GRAVITY [G] <span id="G-val">10</span></div> <!-- Label for Gravitational Constant (G) -->
    <input type="range" id="G" min="0" max="600" step="any" value="10"> <!-- Slider for G: 0 to 600, allows float values -->
    <div class="lbl">SOFTENING [e] <span id="soft-val">1.5</span></div> <!-- Label for softening parameter (prevents singularity division-by-zero errors) -->
    <input type="range" id="soft" min="0" max="25" step="0.1" value="1.5"> <!-- Slider for softening: 0 to 25 -->
    <div class="lbl">TIME STEP [dt] <span id="dt-val">0.002</span></div> <!-- Label for simulation delta time multiplier -->
    <input type="range" id="dt-sl" min="0.0001" max="0.0667" step="any" value="0.002"> <!-- Slider for time step manipulation -->
    <div class="lbl">DRAG COEFF <span id="damp-val">1.00</span></div> <!-- Label for velocity damping (friction/drag) -->
    <input type="range" id="damp" min="0.95" max="1.0" step="0.0001" value="1.0"> <!-- Slider for damping: 1.0 is no drag, <1 is drag -->
    <div class="lbl">MAX VELOCITY <span id="vmax-val">8000</span></div> <!-- Label for terminal velocity cap -->
    <input type="range" id="vmax" min="0" max="10000" step="any" value="8000"> <!-- Slider to prevent bodies from moving infinitely fast -->
  </div> <!-- Close Physics Engine section -->

  <div class="sec"> <!-- Settings group: Relativistic effects -->
    <span class="sec-t">RELATIVITY</span> <!-- Label breaking top border -->
    <label class="chk"><input type="checkbox" id="use-rel" checked> ENABLE RELATIVITY</label> <!-- Toggle for General Relativity equations -->
    <label class="chk"><input type="checkbox" id="use-1pn"> 1PN / EIH CORRECTIONS</label> <!-- Toggle for First Post-Newtonian (Einstein-Infeld-Hoffmann) math -->
    <label class="chk"><input type="checkbox" id="use-2pn"> 2PN HIGHER-ORDER TERMS</label> <!-- Toggle for an extra strong-field correction layered on 1PN -->
    <label class="chk"><input type="checkbox" id="gw-damp" checked> GW DAMPING (HEAVY)</label> <!-- Toggle for simple gravitational-wave energy loss on heavy binaries -->
    <label class="chk"><input type="checkbox" id="show-hor" checked> SHOW EVENT HORIZONS</label> <!-- Toggle for rendering black hole Schwarzschild radii -->
    <div class="lbl">LIGHTSPEED [c] <span id="c-val">8000</span></div> <!-- Label for the speed of light parameter -->
    <input type="range" id="c-sl" min="200" max="10000" step="10" value="8000"> <!-- Slider to adjust the local speed of light -->
    <div class="info" id="rel-info">MODE: PACZYNSKI-WIITA</div> <!-- Note indicating the specific pseudo-Newtonian potential model in use -->
  </div> <!-- Close Relativity section -->

  <div class="sec"> <!-- Settings group: Collision behavior -->
    <span class="sec-t">COLLISION</span> <!-- Label breaking top border -->
    <label class="chk"><input type="checkbox" id="use-frag"> FRAGMENTATION</label> <!-- Toggle allowing bodies to shatter upon high-speed impact -->
    <label class="chk"><input type="checkbox" id="tidal-frag" checked> TIDAL DISRUPTION</label> <!-- Toggle Roche-limit breakup of smaller bodies -->
    <label class="chk"><input type="checkbox" id="ejecta"> MERGER EJECTA</label> <!-- Toggle debris spawned from energetic mergers -->
    <div class="lbl">RESTITUTION [e] <span id="rest-val">0.00</span></div> <!-- Label for Coefficient of Restitution (bounciness) -->
    <input type="range" id="rest-sl" min="0" max="0.99" step="0.01" value="0"> <!-- Slider for elasticity -->
    <div class="info">0=MERGE &nbsp; 1=ELASTIC</div> <!-- Helper text explaining the restitution slider bounds -->
    <div class="lbl">FRAG SPEED <span id="frag-val">250</span></div> <!-- Label for impact velocity required to trigger fragmentation -->
    <input type="range" id="frag-sl" min="10" max="2000" step="10" value="250"> <!-- Slider setting the fragmentation speed threshold -->
  </div> <!-- Close Collision section -->

  <div class="sec"> <!-- Settings group: Optional physics extensions -->
    <span class="sec-t">ADVANCED</span>
    <label class="chk"><input type="checkbox" id="adaptive-dt" checked> ADAPTIVE DT</label>
    <label class="chk"><input type="checkbox" id="use-leapfrog"> LEAPFROG INTEGRATOR</label>
    <label class="chk"><input type="checkbox" id="realistic-radii"> REALISTIC RADII</label>
    <label class="chk"><input type="checkbox" id="dark-matter"> DARK MATTER HALO</label>
  </div>

  <div class="sec"> <!-- Settings group: Execution state controls -->
    <span class="sec-t">REALITY CONTROLS</span> <!-- Label breaking top border -->
    <button id="run-btn">&#9654;&#160; START</button> <!-- Button to initiate the continuous simulation loop -->
    <button id="step-btn">&#9658;&#9658; SINGLE STEP</button> <!-- Button to manually push the physics engine forward 1 tick -->
    <button id="clear-btn" class="warn">&#10005;&#160;&#160; CLEAR </button> <!-- Button to wipe the board; uses .warn class for orange tint -->
    <button id="spawn-btn">&#9675;&#160;&#160; SPAWN 2K</button> <!-- Button to randomly distribute 2000 bodies on the canvas -->
  </div> <!-- Close Reality Controls section -->

  <div class="sec"> <!-- Settings group: Live readout statistics -->
    <span class="sec-t">TELEMETRY</span> <!-- Label breaking top border -->
    <div class="tele" id="stats">BODIES&#160;&#160;0<br>RENDER&#160;&#160;-- fps<br>STEP&#160;&#160;&#160;-- ms<br>MODE&#160;&#160;&#160;NORMAL</div> <!-- Readout lines updated continuously via frontend JS -->
  </div> <!-- Close Telemetry section -->
</div> <!-- Close sidebar controls container -->

<script>
// ── canvas setup ──────────────────────────────────────────────────── // Section: Canvas Initialization
const canvas = document.getElementById("sim");                          // Grab the HTML canvas element by its ID
const ctx    = canvas.getContext("2d");                                 // Get the 2D rendering context for drawing operations
const $      = id => document.getElementById(id);                       // Helper function: jQuery-style DOM element selector

function resizeCanvas() {                                               // Function to dynamically size the canvas to its container
    canvas.width  = canvas.offsetWidth  || 800;                         // Set internal render width to physical CSS width (fallback 800)
    canvas.height = canvas.offsetHeight || 600;                         // Set internal render height to physical CSS height (fallback 600)
}
window.addEventListener("resize", resizeCanvas);                        // Re-calculate dimensions whenever the browser window is resized
resizeCanvas();                                                         // Trigger initial sizing on load

// ── client-side view state ────────────────────────────────────────── // Section: Frontend State Management
const VIEW_SCALE  = 10.0;                                               // Zoom level: how many world units fit in a screen pixel
const THROW_SCALE = 6.0;                                                // Multiplier for mouse-drag velocity when spawning bodies
const HOR_COLOR   = "#cc44ff";                                          // Color for rendering Event Horizons (purple)
const HUD_COLOR   = "#003d33";                                          // Color for the Heads-Up Display text and corners (dark green)
const MASS_MIN    = 1, MASS_MAX = 1000000;                              // Hard limits for particle mass to prevent rendering errors
const R_MIN = 1.0, R_MAX = 10.0;                                        // Visual radius limits (min 1px, max 10px)
const CBRT_MIN = Math.cbrt(MASS_MIN), CBRT_MAX = Math.cbrt(MASS_MAX);   // Precompute cube roots for volume-to-radius scaling

let latestBodies  = [];                                                 // Array holding the most recent physics state from the server
let cameraX = 0, cameraY = 0;                                           // Camera coordinates focusing on the first body
let running = false;                                                    // Flag tracking if the continuous simulation loop is active
let spawnMass = 10000;                                                  // Current selected mass for new bodies
let showHor = true;                                                     // Flag to render black hole event horizons
let useRelLocal = true;                                                 // Mirrored flag for General Relativity to render r_s (Schwarzschild radius) locally
let realisticRadiiLocal = false;                                        // Mirror backend radius mode for local preview rendering

// local param mirrors (for r_s display, horizon display)               // Cache backend physics constants to avoid redundant network calls during rendering
let localG = 10, localC = 8000;                                         // Gravitational constant and Speed of Light

function massToRadius(m) {                                              // Calculate visual size based on mass (assuming constant density, Vol ∝ Mass)
    const cl = Math.min(Math.max(m, MASS_MIN), MASS_MAX);               // Clamp mass within defined bounds
    if (realisticRadiiLocal) {
        const ROCKY_MAX = 100000;
        const GAS_MAX = 1000000;
        if (cl <= ROCKY_MAX) return Math.min(R_MAX, Math.max(R_MIN, Math.cbrt(cl)));
        if (cl <= GAS_MAX) {
            const rockyR = Math.cbrt(ROCKY_MAX);
            return Math.min(R_MAX, Math.max(R_MIN, rockyR * Math.pow(cl / ROCKY_MAX, 0.55)));
        }
        return Math.min(R_MAX, Math.max(R_MIN, 10.0 * Math.cbrt(GAS_MAX / cl)));
    }
    return R_MIN + (Math.cbrt(cl) - CBRT_MIN) * (R_MAX - R_MIN) / (CBRT_MAX - CBRT_MIN); // Linear interpolation on the cube root scale
}
function rsRadius(G, m, c) { return c > 0 ? (2*G*m)/(c*c) : 0; }        // Calculate Schwarzschild radius (event horizon) using local constants
function fmt(n) {                                                       // UI formatting helper for numbers
    const a = Math.abs(n);                                              // Get absolute value
    return a >= 1000 ? n.toFixed(0) : a >= 10 ? n.toFixed(1) : n.toPrecision(3); // Scale decimal precision based on magnitude
}

const DAMPING_REF_DT = 1 / 360;
const MAX_COLL = 8;

function randBodyColor() {
    const ch = () => 140 + Math.floor(Math.random() * 116);
    return `rgb(${ch()},${ch()},${ch()})`;
}

class Sim {
    constructor() {
        this.bodies = [];
        this.next_id = 0;
        this.last_dt = 0.002;
        this.G = 10.0;
        this.softening = 1.5;
        this.dt = 0.002;
        this.damping = 1.0;
        this.max_speed = 8000.0;
        this.c = 8000.0;
        this.use_rel = true;
        this.use_1pn = false;
        this.use_2pn = false;
        this.gw_damping = true;
        this.gw_mass_threshold = 100000.0;
        this.tidal_frag = true;
        this.adaptive_dt = true;
        this.dt_min = 0.0005;
        this.dt_max = 0.01;
        this.realistic_radii = false;
        this.use_leapfrog = false;
        this.ejecta_enabled = false;
        this.dark_matter = false;
        this.dm_scale = 1000.0;
        this.dm_mass = 1e7;
        this.restitution = 0.0;
        this.frag_enabled = false;
        this.frag_threshold = 250.0;
        this.collision_limit = 250;
        this.max_bodies = 1200;
    }

    massToRadius(m) {
        const clamped = Math.min(Math.max(Number(m), MASS_MIN), MASS_MAX);
        if (!this.realistic_radii) {
            return R_MIN + (Math.cbrt(clamped) - CBRT_MIN) * (R_MAX - R_MIN) / (CBRT_MAX - CBRT_MIN);
        }
        const rockyMax = 100000;
        const gasMax = 1000000;
        if (clamped <= rockyMax) return Math.min(R_MAX, Math.max(R_MIN, Math.cbrt(clamped)));
        if (clamped <= gasMax) {
            const rockyR = Math.cbrt(rockyMax);
            return Math.min(R_MAX, Math.max(R_MIN, rockyR * Math.pow(clamped / rockyMax, 0.55)));
        }
        return Math.min(R_MAX, Math.max(R_MIN, 10.0 * Math.cbrt(gasMax / clamped)));
    }

    refreshBodyRadii() {
        for (const b of this.bodies) b.r = this.massToRadius(b.m);
    }

    estimateDt(bs) {
        const dt = Math.max(this.dt_min, Math.min(this.dt_max, this.dt));
        if (!this.adaptive_dt || bs.length < 2) return dt;
        const soft2 = Math.max(this.softening * this.softening, 1e-9);
        const totalMass = bs.reduce((sum, b) => sum + b.m, 0);
        let maxA = 0;
        for (const b of bs) {
            const aEst = this.G * Math.max(totalMass - b.m, 0) / soft2;
            if (aEst > maxA) maxA = aEst;
        }
        if (maxA <= 0) return dt;
        const dtEst = 0.2 * Math.sqrt(soft2 / maxA);
        return Math.max(this.dt_min, Math.min(this.dt_max, Math.min(dt, dtEst)));
    }

    dmAcceleration(x, y) {
        const r2 = x * x + y * y;
        if (r2 < 1e-9) return [0, 0];
        const scale2 = this.dm_scale * this.dm_scale;
        const factor = -this.G * this.dm_mass / Math.pow(r2 + scale2, 1.5);
        return [x * factor, y * factor];
    }

    addBody(x, y, vx, vy, m, color = null, forcedId = null) {
        if (this.bodies.length >= this.max_bodies) return false;
        const id = forcedId == null ? this.next_id : forcedId;
        this.next_id = Math.max(this.next_id, id + 1);
        this.bodies.push({
            id,
            x,
            y,
            vx,
            vy,
            m: Number(m),
            r: this.massToRadius(m),
            color: color || randBodyColor(),
        });
        return true;
    }

    loadBodies(bodies) {
        this.bodies = [];
        this.next_id = 0;
        for (const b of bodies || []) {
            this.addBody(
                Number(b.x),
                Number(b.y),
                Number(b.vx),
                Number(b.vy),
                Number(b.m),
                b.color || randBodyColor(),
                Number.isFinite(b.id) ? Number(b.id) : null,
            );
        }
    }

    loadSnapshot(state) {
        if (!state) return;
        if (state.params) this.setParams(state.params);
        if (Array.isArray(state.bodies)) this.loadBodies(state.bodies);
    }

    tidalFragment(i, j) {
        let bi = this.bodies[i];
        let bj = this.bodies[j];
        if (bi.m > bj.m) {
            [i, j] = [j, i];
            bi = this.bodies[i];
            bj = this.bodies[j];
        }
        if (!bi || !bj || bi.m <= 0) return false;

        const roche = bi.r * Math.pow((2 * bj.m) / bi.m, 1 / 3);
        const dx = bj.x - bi.x;
        const dy = bj.y - bi.y;
        const dist = Math.hypot(dx, dy);
        if (dist <= 0 || dist >= roche) return false;

        const disrupted = this.bodies.splice(i, 1)[0];
        const fragMass = Math.max(1.0, disrupted.m / 3.0);
        const fr = this.massToRadius(fragMass);
        const tangent = Math.atan2(dy, dx) + Math.PI / 2;
        const relSpeed = Math.hypot(disrupted.vx - bj.vx, disrupted.vy - bj.vy);
        const spread = Math.max(relSpeed, this.frag_threshold * 0.15);
        for (let k = 0; k < 3; k++) {
            const angle = tangent + (k - 1) * 0.65 + (Math.random() - 0.5) * 0.2;
            const speed = spread * (0.6 + Math.random() * 0.4);
            this.addBody(
                disrupted.x + Math.cos(angle) * (fr + 1) * 2.0,
                disrupted.y + Math.sin(angle) * (fr + 1) * 2.0,
                disrupted.vx + Math.cos(angle) * speed,
                disrupted.vy + Math.sin(angle) * speed,
                fragMass,
            );
        }
        return true;
    }

    merge(si, di) {
        const s = this.bodies[si];
        const d = this.bodies[di];
        const m1 = s.m;
        const m2 = d.m;
        const mt = m1 + m2 || 1;
        const vrelx = s.vx - d.vx;
        const vrely = s.vy - d.vy;
        const relKe = 0.5 * (m1 * m2 / mt) * (vrelx * vrelx + vrely * vrely);
        let ejectaMass = 0;
        if (this.ejecta_enabled && this.restitution < 0.05 && mt > 2.0 && relKe > this.frag_threshold * 0.5) {
            ejectaMass = Math.min(Math.max(1.0, mt * 0.01), mt * 0.25);
        }
        const mergedMass = mt - ejectaMass;
        d.x = (m1 * s.x + m2 * d.x) / mt;
        d.y = (m1 * s.y + m2 * d.y) / mt;
        d.vx = (m1 * s.vx + m2 * d.vx) / mt;
        d.vy = (m1 * s.vy + m2 * d.vy) / mt;
        d.m = mergedMass;
        d.r = this.massToRadius(mergedMass);
        this.bodies.splice(si, 1);
        if (ejectaMass > 0) {
            const angle = Math.random() * Math.PI * 2;
            const speed = Math.sqrt(2 * relKe / ejectaMass) * 0.25;
            this.addBody(
                d.x,
                d.y,
                d.vx + Math.cos(angle) * speed,
                d.vy + Math.sin(angle) * speed,
                ejectaMass,
            );
        }
    }

    fragment(si, di) {
        const s = this.bodies[si];
        const d = this.bodies[di];
        const m1 = s.m;
        const m2 = d.m;
        const mt = m1 + m2;
        const cmx = (m1 * s.x + m2 * d.x) / mt;
        const cmy = (m1 * s.y + m2 * d.y) / mt;
        const cmvx = (m1 * s.vx + m2 * d.vx) / mt;
        const cmvy = (m1 * s.vy + m2 * d.vy) / mt;
        const dvx = s.vx - d.vx;
        const dvy = s.vy - d.vy;
        const relKe = 0.5 * (m1 * m2 / mt) * (dvx * dvx + dvy * dvy);
        const ejSpeed = Math.sqrt(Math.max(0, 2 * relKe / mt)) * 0.5 + this.frag_threshold * 0.15;
        const hi = Math.max(si, di);
        const lo = Math.min(si, di);
        this.bodies.splice(hi, 1);
        this.bodies.splice(lo, 1);
        const fragMass = Math.max(1.0, mt / 3);
        const fr = this.massToRadius(fragMass);
        for (let k = 0; k < 3; k++) {
            const angle = 2 * Math.PI * k / 3 + (Math.random() - 0.5) * 0.4;
            const speed = ejSpeed * (0.8 + Math.random() * 0.4);
            this.addBody(
                cmx + Math.cos(angle) * (fr + 1) * 2.5,
                cmy + Math.sin(angle) * (fr + 1) * 2.5,
                cmvx + Math.cos(angle) * speed,
                cmvy + Math.sin(angle) * speed,
                fragMass,
            );
        }
    }

    collectPairs() {
        const out = [];
        for (let i = 0; i < this.bodies.length - 1; i++) {
            const bi = this.bodies[i];
            for (let j = i + 1; j < this.bodies.length; j++) {
                const bj = this.bodies[j];
                const dx = bj.x - bi.x;
                const dy = bj.y - bi.y;
                const lim = bi.r + bj.r;
                const d2 = dx * dx + dy * dy;
                if (d2 >= lim * lim) continue;
                const dist = Math.sqrt(d2);
                out.push([lim - dist, bi.id * 2000000 + bj.id, i, j, dist, dx, dy]);
            }
        }
        out.sort((a, b) => (b[0] - a[0]) || (b[1] - a[1]));
        return out;
    }

    resolve() {
        for (let iter = 0; iter < MAX_COLL; iter++) {
            const pairs = this.collectPairs();
            if (!pairs.length) break;
            const [, , i, j, dist, dx, dy] = pairs[0];
            const bi = this.bodies[i];
            const bj = this.bodies[j];
            const d = Math.max(dist, 1e-9);
            const nx = dx / d;
            const ny = dy / d;
            const dvx = bi.vx - bj.vx;
            const dvy = bi.vy - bj.vy;
            const dvn = dvx * nx + dvy * ny;

            if (this.tidal_frag && this.restitution < 0.05 && this.tidalFragment(i, j)) {
                continue;
            }
            if (this.restitution < 0.05) {
                if (this.frag_enabled && Math.abs(dvn) > this.frag_threshold) this.fragment(i, j);
                else {
                    const dst = bi.m >= bj.m ? i : j;
                    const src = dst === i ? j : i;
                    this.merge(src, dst);
                }
                continue;
            }
            if (dvn > 0) {
                const mi = bi.m;
                const mj = bj.m;
                const ji = (1 + this.restitution) * dvn * mi * mj / (mi + mj);
                bi.vx -= ji / mi * nx;
                bi.vy -= ji / mi * ny;
                bj.vx += ji / mj * nx;
                bj.vy += ji / mj * ny;
            }
            const sep = bi.r + bj.r - dist + 0.05;
            const mt = bi.m + bj.m;
            bi.x -= nx * sep * bj.m / mt;
            bi.y -= ny * sep * bj.m / mt;
            bj.x += nx * sep * bi.m / mt;
            bj.y += ny * sep * bi.m / mt;
        }
    }

    computeAccelerations(bs, dt) {
        const ax = new Array(bs.length).fill(0);
        const ay = new Array(bs.length).fill(0);
        const G = this.G;
        const c = this.c;
        const eps2 = this.softening * this.softening;

        for (let i = 0; i < bs.length - 1; i++) {
            const bi = bs[i];
            for (let j = i + 1; j < bs.length; j++) {
                const bj = bs[j];
                const dx = bj.x - bi.x;
                const dy = bj.y - bi.y;
                const rs2 = dx * dx + dy * dy + eps2;
                if (rs2 < 1e-9) continue;

                if (!this.use_rel) {
                    const ir3 = Math.pow(rs2, -1.5);
                    const si = G * bj.m * ir3;
                    const sj = G * bi.m * ir3;
                    ax[i] += dx * si;
                    ay[i] += dy * si;
                    ax[j] -= dx * sj;
                    ay[j] -= dy * sj;
                } else if (!this.use_1pn) {
                    const rs = Math.sqrt(rs2);
                    const di = Math.max(rs - rsRadius(G, bj.m, c), 1);
                    const dj = Math.max(rs - rsRadius(G, bi.m, c), 1);
                    const ir = 1 / rs;
                    const si = G * bj.m * ir / (di * di);
                    const sj = G * bi.m * ir / (dj * dj);
                    ax[i] += dx * si;
                    ay[i] += dy * si;
                    ax[j] -= dx * sj;
                    ay[j] -= dy * sj;
                } else {
                    const rs = Math.sqrt(rs2);
                    const ir = 1 / rs;
                    const ir3 = ir / rs2;
                    const si = G * bj.m * ir3;
                    const sj = G * bi.m * ir3;
                    ax[i] += dx * si;
                    ay[i] += dy * si;
                    ax[j] -= dx * sj;
                    ay[j] -= dy * sj;

                    const c2 = c * c;
                    const nx = dx * ir;
                    const ny = dy * ir;
                    const vi2 = bi.vx * bi.vx + bi.vy * bi.vy;
                    const vj2 = bj.vx * bj.vx + bj.vy * bj.vy;
                    const vivj = bi.vx * bj.vx + bi.vy * bj.vy;
                    const vi_n = bi.vx * nx + bi.vy * ny;
                    const vj_n = bj.vx * nx + bj.vy * ny;
                    const Gmi_r = G * bi.m * ir;
                    const Gmj_r = G * bj.m * ir;
                    const ir2 = ir * ir;
                    const dvx2 = bi.vx - bj.vx;
                    const dvy2 = bi.vy - bj.vy;

                    const fi_f = G * bj.m * ir2 / c2;
                    const fi_s = -vi2 - 2 * vj2 + 4 * vivj + 1.5 * vj_n * vj_n + 5 * Gmi_r + 4 * Gmj_r;
                    const fi_v = 4 * vi_n - 3 * vj_n;
                    const pn1_ix = fi_f * (nx * fi_s + dvx2 * fi_v);
                    const pn1_iy = fi_f * (ny * fi_s + dvy2 * fi_v);
                    ax[i] += pn1_ix;
                    ay[i] += pn1_iy;

                    const fj_f = G * bi.m * ir2 / c2;
                    const fj_s = -vj2 - 2 * vi2 + 4 * vivj + 1.5 * vi_n * vi_n + 5 * Gmj_r + 4 * Gmi_r;
                    const fj_v = 4 * vj_n - 3 * vi_n;
                    const pn1_jx = fj_f * (-nx * fj_s + dvx2 * fj_v);
                    const pn1_jy = fj_f * (-ny * fj_s + dvy2 * fj_v);
                    ax[j] += pn1_jx;
                    ay[j] += pn1_jy;

                    if (this.use_2pn) {
                        const pnBoost = Math.min(0.35, 2.0 * (0.5 * (vi2 + vj2) + Gmi_r + Gmj_r) / c2);
                        ax[i] += pnBoost * pn1_ix;
                        ay[i] += pnBoost * pn1_iy;
                        ax[j] += pnBoost * pn1_jx;
                        ay[j] += pnBoost * pn1_jy;
                    }
                }

                if (this.use_rel && this.gw_damping && (bi.m >= this.gw_mass_threshold || bj.m >= this.gw_mass_threshold)) {
                    const r = Math.sqrt(rs2);
                    if (r > 0) {
                        const vrelx = bj.vx - bi.vx;
                        const vrely = bj.vy - bi.vy;
                        const vrel2 = vrelx * vrelx + vrely * vrely;
                        if (vrel2 > 1e-12) {
                            const vrel = Math.sqrt(vrel2);
                            const ux = vrelx / vrel;
                            const uy = vrely / vrel;
                            const beta2 = Math.min(1.0, vrel2 / Math.max(c * c, 1e-9));
                            const mi = Math.max(bi.m, 1e-9);
                            const mj = Math.max(bj.m, 1e-9);
                            const dE_dt = (32.0 / 5.0) * (
                                Math.pow(G, 4) * mi * mi * mj * mj * (mi + mj) / Math.max(Math.pow(c, 5) * Math.pow(r, 5), 1e-9)
                            );
                            let force = dE_dt / vrel;
                            force *= beta2;
                            const reducedMass = mi * mj / (mi + mj);
                            const maxForce = 0.25 * reducedMass * vrel / Math.max(dt, 1e-9);
                            force = Math.min(force, maxForce);
                            ax[i] += ux * (force / mi);
                            ay[i] += uy * (force / mi);
                            ax[j] -= ux * (force / mj);
                            ay[j] -= uy * (force / mj);
                        }
                    }
                }
            }
        }

        if (this.dark_matter) {
            for (let i = 0; i < bs.length; i++) {
                const [dax, day] = this.dmAcceleration(bs[i].x, bs[i].y);
                ax[i] += dax;
                ay[i] += day;
            }
        }
        return {ax, ay};
    }

    stepOnce() {
        if (!this.bodies.length) return;
        const doCollisions = this.bodies.length <= this.collision_limit;
        if (doCollisions) this.resolve();
        const dt = this.estimateDt(this.bodies);
        this.last_dt = dt;
        const drag = Math.pow(this.damping, dt / DAMPING_REF_DT);
        const vcap = this.use_rel ? Math.min(this.max_speed, this.c) : this.max_speed;
        const vcap2 = vcap * vcap;

        if (this.use_leapfrog) {
            let {ax, ay} = this.computeAccelerations(this.bodies, dt);
            for (let i = 0; i < this.bodies.length; i++) {
                this.bodies[i].vx += ax[i] * dt * 0.5;
                this.bodies[i].vy += ay[i] * dt * 0.5;
            }
            for (const b of this.bodies) {
                b.x += b.vx * dt;
                b.y += b.vy * dt;
            }
            ({ax, ay} = this.computeAccelerations(this.bodies, dt));
            for (let i = 0; i < this.bodies.length; i++) {
                const b = this.bodies[i];
                b.vx = (b.vx + ax[i] * dt * 0.5) * drag;
                b.vy = (b.vy + ay[i] * dt * 0.5) * drag;
                const v2 = b.vx * b.vx + b.vy * b.vy;
                if (v2 > vcap2) {
                    const scale = vcap / Math.sqrt(v2);
                    b.vx *= scale;
                    b.vy *= scale;
                }
            }
        } else {
            const {ax, ay} = this.computeAccelerations(this.bodies, dt);
            for (let i = 0; i < this.bodies.length; i++) {
                const b = this.bodies[i];
                b.vx = (b.vx + ax[i] * dt) * drag;
                b.vy = (b.vy + ay[i] * dt) * drag;
                const v2 = b.vx * b.vx + b.vy * b.vy;
                if (v2 > vcap2) {
                    const scale = vcap / Math.sqrt(v2);
                    b.vx *= scale;
                    b.vy *= scale;
                }
                b.x += b.vx * dt;
                b.y += b.vy * dt;
            }
        }
        if (doCollisions) this.resolve();
    }

    step(n = 1) {
        for (let i = 0; i < n; i++) this.stepOnce();
    }

    setParams(params) {
        const floatKeys = new Set(["G", "softening", "dt", "damping", "max_speed", "c", "frag_threshold", "restitution", "dt_min", "dt_max", "dm_scale", "dm_mass"]);
        const boolKeys = new Set(["use_rel", "use_1pn", "use_2pn", "gw_damping", "tidal_frag", "adaptive_dt", "realistic_radii", "use_leapfrog", "ejecta_enabled", "dark_matter", "frag_enabled"]);
        let refreshRadii = false;
        for (const [key, value] of Object.entries(params || {})) {
            if (floatKeys.has(key)) this[key] = Number(value);
            else if (boolKeys.has(key)) this[key] = !!value;
            if (key === "realistic_radii") refreshRadii = true;
        }
        if (this.dt_min > this.dt_max) [this.dt_min, this.dt_max] = [this.dt_max, this.dt_min];
        if (this.use_2pn) this.use_1pn = true;
        else if (!this.use_1pn) this.use_2pn = false;
        if (refreshRadii) this.refreshBodyRadii();
    }

    params() {
        return {
            G: this.G,
            softening: this.softening,
            dt: this.dt,
            damping: this.damping,
            max_speed: this.max_speed,
            c: this.c,
            use_rel: this.use_rel,
            use_1pn: this.use_1pn,
            use_2pn: this.use_2pn,
            gw_damping: this.gw_damping,
            tidal_frag: this.tidal_frag,
            adaptive_dt: this.adaptive_dt,
            dt_min: this.dt_min,
            dt_max: this.dt_max,
            realistic_radii: this.realistic_radii,
            use_leapfrog: this.use_leapfrog,
            ejecta_enabled: this.ejecta_enabled,
            dark_matter: this.dark_matter,
            dm_scale: this.dm_scale,
            dm_mass: this.dm_mass,
            restitution: this.restitution,
            frag_enabled: this.frag_enabled,
            frag_threshold: this.frag_threshold,
        };
    }

    clear() {
        this.bodies = [];
    }
}

const sim = new Sim();
latestBodies = sim.bodies;

function syncViewFromSim() {
    latestBodies = sim.bodies;
    if (latestBodies.length > 0) {
        cameraX = latestBodies[0].x;
        cameraY = latestBodies[0].y;
    }
}

function syncUiFromSim() {
    const params = sim.params();
    localG = params.G;
    localC = params.c;
    useRelLocal = !!params.use_rel;
    realisticRadiiLocal = !!params.realistic_radii;
    $("G").value = params.G; $("G-val").textContent = fmt(params.G);
    $("soft").value = params.softening; $("soft-val").textContent = fmt(params.softening);
    $("dt-sl").value = params.dt; $("dt-val").textContent = fmt(params.dt);
    $("damp").value = params.damping; $("damp-val").textContent = fmt(params.damping);
    $("vmax").value = params.max_speed; $("vmax-val").textContent = fmt(params.max_speed);
    $("c-sl").value = params.c; $("c-val").textContent = fmt(params.c);
    $("rest-sl").value = params.restitution; $("rest-val").textContent = Number(params.restitution).toFixed(2);
    $("frag-sl").value = params.frag_threshold; $("frag-val").textContent = fmt(params.frag_threshold);
    $("use-rel").checked = !!params.use_rel;
    $("use-1pn").checked = !!params.use_1pn;
    $("use-2pn").checked = !!params.use_2pn;
    $("gw-damp").checked = !!params.gw_damping;
    $("use-frag").checked = !!params.frag_enabled;
    $("tidal-frag").checked = !!params.tidal_frag;
    $("ejecta").checked = !!params.ejecta_enabled;
    $("adaptive-dt").checked = !!params.adaptive_dt;
    $("use-leapfrog").checked = !!params.use_leapfrog;
    $("realistic-radii").checked = !!params.realistic_radii;
    $("dark-matter").checked = !!params.dark_matter;
    updateMassInfo();
    updateRelInfo();
}

async function fetchInitialState() {
    try {
        const r = await fetch('/state');
        if (!r.ok) return;
        applyState(await r.json());
    } catch (err) {
        console.warn('Initial state load failed; starting locally.', err);
    }
}

// ── API helpers ───────────────────────────────────────────────────── // Section: Backend Communication layer
function applyState(state) {
    if (!state) return;
    sim.loadSnapshot(state);
    syncUiFromSim();
    syncViewFromSim();
}

function applyParams(params) {
    if (!params) return;
    sim.setParams(params);
    syncUiFromSim();
    syncViewFromSim();
}

function appendSpawnedBodies(bodies) {                                  // Locally append bodies before the server responds (optimistic UI update)
    let added = 0;
    if (!Array.isArray(bodies) || bodies.length === 0) return added;
    for (const b of bodies) {
        if (!sim.addBody(b.x, b.y, b.vx, b.vy, b.m)) break;
        added += 1;
    }
    syncViewFromSim();
    return added;
}

function nextFrame() {                                                  // Helper to yield execution to the browser's render pipeline
    return new Promise(resolve => requestAnimationFrame(resolve));      // Returns a promise that resolves on the next frame (prevents UI freezing during loops)
}

// ── step loop (decoupled from render) ─────────────────────────────── // Section: Physics Tick Management (Asynchronous to framerate)
let stepMs = null;

function stepsPerFrame(bodyCount) {
    if (bodyCount <= 80) return 4;
    if (bodyCount <= 150) return 2;
    if (bodyCount <= 300) return 1;
    return 0;
}

function runPhysicsFrame() {
    if (!running) return;
    const steps = stepsPerFrame(latestBodies.length);
    if (steps <= 0) {
        stepMs = 0;
        return;
    }
    const t0 = performance.now();
    sim.step(steps);
    stepMs = performance.now() - t0;
    syncViewFromSim();
}

// ── render loop ───────────────────────────────────────────────────── // Section: Main Render Pipeline
let lastT = performance.now(), fpsVal = null, frameCount = 0, fpsAccum = 0; // Telemetry tracking variables

function worldToScreen(x, y) {                                          // Transform world coordinates to screen pixel coordinates
    return {
        x: (x - cameraX) / VIEW_SCALE + canvas.width  / 2,              // Offset by camera, scale down, and center on screen X
        y: (y - cameraY) / VIEW_SCALE + canvas.height / 2,              // Offset by camera, scale down, and center on screen Y
    };
}
function screenToWorld(x, y) {                                          // Transform screen pixel coordinates back to world coordinates
    return {
        x: (x - canvas.width  / 2) * VIEW_SCALE + cameraX,              // Un-center, scale up, add camera offset X
        y: (y - canvas.height / 2) * VIEW_SCALE + cameraY,              // Un-center, scale up, add camera offset Y
    };
}

function drawHUD() {                                                    // Render the Heads-Up Display
    const W = canvas.width, H = canvas.height, L = 30;                  // Setup dimensions and corner line length
    ctx.strokeStyle = HUD_COLOR; ctx.lineWidth = 2; ctx.setLineDash([]); // Set stroke style for corners
    for (const [ox,oy,sx,sy] of [[4,4,1,1],[W-4,4,-1,1],[4,H-4,1,-1],[W-4,H-4,-1,-1]]) { // Iterate over the 4 corners
        ctx.beginPath();
        ctx.moveTo(ox,oy); ctx.lineTo(ox+sx*L,oy);                      // Draw horizontal corner segment
        ctx.moveTo(ox,oy); ctx.lineTo(ox,oy+sy*L);                      // Draw vertical corner segment
        ctx.stroke();                                                   // Execute stroke
    }
    ctx.fillStyle = HUD_COLOR; ctx.font = "11px 'Courier New'";         // Setup text style
    ctx.fillText("STELLAR DYNAMICS  //  GRAVITY ENGINE  [JS]", 16, 20); // Render title text
    ctx.setLineDash([]);                                                // Reset line dash for subsequent renders
}

function render(now) {                                                  // Main Animation Frame callback
    const elapsed = Math.min((now - lastT) / 1000, 0.25);               // Calculate delta-time (cap at 250ms to prevent huge jumps if tab is backgrounded)
    lastT = now;                                                        // Store timestamp for next frame
    frameCount++; fpsAccum += elapsed;                                  // Increment telemetry counters
    if (fpsAccum >= 0.5) {                                              // Every 0.5 seconds...
        fpsVal = frameCount / fpsAccum;                                 // Calculate average FPS
        frameCount = fpsAccum = 0;                                      // Reset counters
        updateStats();                                                  // Push new telemetry to the UI DOM
    }

    runPhysicsFrame();

    ctx.fillStyle = "#050a10";                                          // Set deep space background color
    ctx.fillRect(0, 0, canvas.width, canvas.height);                    // Clear previous frame by filling the whole canvas

    if (showHor && useRelLocal) {                                       // If General Relativity and Event Horizons are enabled...
        ctx.strokeStyle = HOR_COLOR; ctx.lineWidth = 1; ctx.setLineDash([3,3]); // Set purple dashed line style
        for (const b of latestBodies) {                                 // Loop through all bodies
            const rs = rsRadius(localG, b.m, localC);                   // Calculate Schwarzschild radius
            if (rs <= 1) continue;                                      // Optimization: don't draw if it's too small to see
            const p = worldToScreen(b.x, b.y);                          // Get screen coordinates
            if (p.x+rs<0||p.x-rs>canvas.width||p.y+rs<0||p.y-rs>canvas.height) continue; // View frustum culling: Skip if off-screen
            ctx.beginPath(); ctx.arc(p.x, p.y, rs, 0, Math.PI*2); ctx.stroke(); // Draw the horizon circle
        }
        ctx.setLineDash([]);                                            // Reset line style
    }

    for (const b of latestBodies) {                                     // Pass 2: Draw the actual celestial bodies
        const p = worldToScreen(b.x, b.y);                              // Get screen coordinates
        if (p.x+b.r<0||p.x-b.r>canvas.width||p.y+b.r<0||p.y-b.r>canvas.height) continue; // View frustum culling: Skip if off-screen
        ctx.fillStyle = b.color;                                        // Use the assigned random color
        ctx.beginPath(); ctx.arc(p.x, p.y, b.r, 0, Math.PI*2); ctx.fill(); // Draw a solid circle
    }

    if (dragging) {                                                     // If the user is currently dragging to spawn a body...
        const start = worldToScreen(dragStartWorld.x, dragStartWorld.y); // Get screen coords of the initial click
        const r = sim.massToRadius(spawnMass);                          // Get the visual size of the pending body
        ctx.strokeStyle = "#00e5cc"; ctx.lineWidth = 1; ctx.setLineDash([]); // Set solid cyan line
        ctx.beginPath(); ctx.arc(start.x, start.y, r, 0, Math.PI*2); ctx.stroke(); // Draw the ghost body where they clicked
        ctx.setLineDash([4,2]);                                         // Set dashed line style for the velocity vector
        ctx.beginPath(); ctx.moveTo(start.x, start.y);                  // Start vector at body origin
        ctx.lineTo(dragCurScreen.x, dragCurScreen.y); ctx.stroke(); ctx.setLineDash([]); // Draw line to current mouse position (represents throw velocity)
    }

    drawHUD();                                                          // Overlay the HUD on top of everything
    requestAnimationFrame(render);                                      // Queue up the next render frame
}
requestAnimationFrame(render);                                          // Kickstart the render loop

// ── mouse ──────────────────────────────────────────────────────────── // Section: Input Handling
let dragging = false, dragStartWorld = {x:0,y:0}, dragCurScreen = {x:0,y:0}; // State variables for drag-to-spawn interactions
function canvasXY(e) {                                                  // Helper to normalize mouse coordinates relative to the canvas element
    const r = canvas.getBoundingClientRect();                           // Get DOM element bounds
    return {
        x: (e.clientX-r.left)*(canvas.width/r.width),                   // Calculate local X accounting for CSS scaling
        y: (e.clientY-r.top)*(canvas.height/r.height),                  // Calculate local Y accounting for CSS scaling
    };
}
canvas.addEventListener("mousedown", e => {                             // On left click down...
    dragCurScreen  = canvasXY(e);                                       // Record current screen coords
    dragStartWorld = screenToWorld(dragCurScreen.x, dragCurScreen.y);   // Lock in the starting world coordinates for the spawn
    dragging = true;                                                    // Flag interaction as active
});
window.addEventListener("mousemove", e => { if (dragging) dragCurScreen = canvasXY(e); }); // Update target vector as mouse moves
window.addEventListener("mouseup", e => {                               // On click release...
    if (!dragging) return;                                              // Ignore if we weren't dragging
    dragging = false;                                                   // End interaction
    dragCurScreen = canvasXY(e);                                        // Get final mouse position
    const p = screenToWorld(dragCurScreen.x, dragCurScreen.y);          // Convert to world coordinates
    sim.addBody(
        dragStartWorld.x,
        dragStartWorld.y,
        (dragStartWorld.x - p.x) * THROW_SCALE,
        (dragStartWorld.y - p.y) * THROW_SCALE,
        spawnMass
    );
    syncViewFromSim();
    updateStats();                                                      // Update body count UI
});
window.addEventListener("blur", () => { dragging = false; });           // Cancel drag if window loses focus

// ── UI ─────────────────────────────────────────────────────────────── // Section: HTML UI Data Binding
function updateMassInfo() {                                             // Update the dynamically calculated UI text for the body spawner
    const r  = sim.massToRadius(spawnMass);                             // Calculate render size
    const rs = useRelLocal ? rsRadius(localG, spawnMass, localC) : 0;   // Calculate Event Horizon size
    $("mass-val").textContent = spawnMass.toLocaleString();             // Update slider label
    $("body-info").innerHTML  =                                         // Inject rich HTML text
        `MASS&#160;&#160;&#160;${spawnMass.toLocaleString()}<br>` +     // Display mass
        `RADIUS&#160;${r.toFixed(2)} px<br>` +                          // Display render radius
        `r_s&#160;&#160;&#160;&#160;${rs.toFixed(2)} px`;               // Display Event horizon radius
}
function updateStats() {                                                // Update the Telemetry panel
    const fps  = fpsVal  != null ? fpsVal.toFixed(1)  : "--";           // Format Frames Per Second
    const sms  = stepMs  != null ? stepMs.toFixed(1)  : "--";           // Format Network+Compute Step latency
    const mode = latestBodies.length >= sim.max_bodies ? "BODY CAP" : (stepsPerFrame(latestBodies.length) === 0 ? "THROTTLED" : "NORMAL");
    $("stats").innerHTML =
        `BODIES&#160;&#160;${latestBodies.length.toLocaleString()}<br>` + // Show count
        `RENDER&#160;&#160;${fps} fps<br>` +
        `STEP&#160;&#160;&#160;${sms} ms<br>` +
        `MODE&#160;&#160;&#160;${mode}`;
}
function updateRelInfo() {                                              // Update the subtext under the Relativity panel
    const rel  = sim.use_rel;
    const pn1  = sim.use_1pn;
    const pn2  = sim.use_2pn;
    const gw   = sim.gw_damping;
    let mode = !rel ? "NEWTONIAN" : pn2 ? "2PN / EIH" : pn1 ? "1PN / EIH" : "PACZYNSKI-WIITA"; // Determine which specific equations the backend is using
    if (rel && gw) mode += " + GW LOSS";
    $("rel-info").textContent = "MODE: " + mode;                        // Update text
}

// param → server helper
function bindParam(id, paramKey, transform, onUpdate) {                 // Utility function to bind an HTML slider directly to the local simulation
    $(id).addEventListener("input", function() {                        // Listen for slider movement
        const v = transform(this.value);                                // Convert HTML string value to integer/float
        applyParams({[paramKey]: v});
        if (onUpdate) onUpdate(v);                                      // Trigger any local UI cleanup callbacks
    });
}

// Bind all HTML sliders to their respective Python variables
bindParam("G",       "G",            parseFloat, v => { localG=v; $("G-val").textContent=fmt(v); updateMassInfo(); });
bindParam("soft",    "softening",    parseFloat, v => $("soft-val").textContent=fmt(v));
bindParam("dt-sl",   "dt",           parseFloat, v => $("dt-val").textContent=fmt(v));
bindParam("damp",    "damping",      parseFloat, v => $("damp-val").textContent=fmt(v));
bindParam("vmax",    "max_speed",    parseFloat, v => $("vmax-val").textContent=fmt(v));
bindParam("c-sl",    "c",            parseFloat, v => { localC=v; $("c-val").textContent=fmt(v); updateMassInfo(); });
bindParam("rest-sl", "restitution",  parseFloat, v => $("rest-val").textContent=v.toFixed(2));
bindParam("frag-sl", "frag_threshold", parseFloat, v => $("frag-val").textContent=fmt(v));

// Setup manual event listeners for items that don't fit the generic bindParam helper
$("mass").addEventListener("input", function() {                        // Spawning mass is purely local until a body is actually spawned
    spawnMass = Math.round(parseFloat(this.value));
    $("mass-val").textContent = spawnMass.toLocaleString();
    updateMassInfo();
});
$("use-rel").addEventListener("change", function() {
    applyParams({use_rel: this.checked});
    updateRelInfo(); updateMassInfo();
});
$("use-1pn").addEventListener("change", function() {
    if (!this.checked && $("use-2pn").checked) $("use-2pn").checked = false;
    applyParams({use_1pn: this.checked, use_2pn: $("use-2pn").checked});
    updateRelInfo();
});
$("use-2pn").addEventListener("change", function() {
    if (this.checked && !$("use-1pn").checked) $("use-1pn").checked = true;
    applyParams({use_1pn: $("use-1pn").checked, use_2pn: this.checked});
    updateRelInfo();
});
$("gw-damp").addEventListener("change", function() {
    applyParams({gw_damping: this.checked});
    updateRelInfo();
});
$("show-hor").addEventListener("change", function() { showHor = this.checked; }); // Local toggle for UI only
$("use-frag").addEventListener("change", function() { applyParams({frag_enabled: this.checked}); });
$("tidal-frag").addEventListener("change", function() { applyParams({tidal_frag: this.checked}); });
$("ejecta").addEventListener("change", function() { applyParams({ejecta_enabled: this.checked}); });
$("adaptive-dt").addEventListener("change", function() { applyParams({adaptive_dt: this.checked}); });
$("use-leapfrog").addEventListener("change", function() { applyParams({use_leapfrog: this.checked}); });
$("realistic-radii").addEventListener("change", function() { applyParams({realistic_radii: this.checked}); updateMassInfo(); });
$("dark-matter").addEventListener("change", function() { applyParams({dark_matter: this.checked}); });

const runBtn = $("run-btn");
function setRunning(nextRunning) {                                      // Manage the master Engine Start/Stop state
    running = nextRunning;
    runBtn.textContent = running ? "|| HOLD ENGINES" : "\u25B6\u00A0 ENGAGE ENGINES"; // Update button text visually
    runBtn.classList.toggle("active", running);                         // Toggle orange highlight class
}

runBtn.addEventListener("click", () => {
    setRunning(!running);                                               // Toggle on click
});

$("step-btn").addEventListener("click", () => {                         // Manual single-step forward
    sim.step(1);
    syncViewFromSim();
    updateStats();
});
$("clear-btn").addEventListener("click", () => {                        // Nuke the board
    sim.clear();
    syncViewFromSim();
    updateStats();
});
$("spawn-btn").addEventListener("click", async () => {                  // Heavy load testing button
    const wasRunning = running;
    setRunning(false);

    const CHUNK_SIZE = 100;
    const TOTAL = 2000;
    const pad   = 50 * VIEW_SCALE;                                      // World padding
    const spanW = canvas.width  * VIEW_SCALE;                           // Calculate visible world width
    const spanH = canvas.height * VIEW_SCALE;                           // Calculate visible world height
    const left  = cameraX - spanW/2 + pad;                              // Find top left world corner
    const top   = cameraY - spanH/2 + pad;
    const spawnBtn = $("spawn-btn");

    spawnBtn.disabled = true;
    try {
        for (let i = 0; i < TOTAL; i += CHUNK_SIZE) {
            let hitCap = false;
            const end = Math.min(i + CHUNK_SIZE, TOTAL);
            for (let j = i; j < end; j++) {
                if (!sim.addBody(
                    left + Math.random() * (spanW - 2 * pad),
                    top  + Math.random() * (spanH - 2 * pad),
                    Math.random() - 0.5,
                    Math.random() - 0.5,
                    spawnMass
                )) {
                    hitCap = true;
                    break;
                }
            }
            syncViewFromSim();
            updateStats();
            if (hitCap) break;
            await nextFrame();
        }
    } finally {
        spawnBtn.disabled = false;
        if (wasRunning) setRunning(true);
    }
});

// init display
$("G-val").textContent    = fmt(10);                                    // Hardcode initial UI values to match backend defaults
$("soft-val").textContent = fmt(1.5);
$("dt-val").textContent   = fmt(0.002);
$("damp-val").textContent = (1.0).toFixed(2);
$("vmax-val").textContent = fmt(8000);
$("c-val").textContent    = fmt(8000);
$("rest-val").textContent = "0.00";
$("frag-val").textContent = fmt(250);
updateMassInfo(); updateRelInfo(); updateStats();                       // Trigger initial UI text generation

// fetch initial state
syncUiFromSim();
syncViewFromSim();
updateStats();
fetchInitialState().then(() => updateStats());
</script>
</body>
</html>"""


# ═══════════════════════════════════════════════════════════════════ # Top visual divider for the Server execution section
# Server                                                              # Section: Standard Python Server Init
# ═══════════════════════════════════════════════════════════════════ # Bottom visual divider
class _Server(HTTPServer):
    allow_reuse_address = True

def main():
    try:
        srv = _Server(("127.0.0.1", PORT), _H)
    except OSError as exc:
        print(f"\n  Failed to start on 127.0.0.1:{PORT}: {exc}\n")
        return 1
    url = f"http://127.0.0.1:{PORT}"
    print(f"\n  Gravity Sim  >  {url}")
    print("  Press Ctrl+C to stop.\n")
    threading.Thread(target=webbrowser.open, args=(url,), daemon=True).start() # Spawn a background daemon thread to automatically open the default web browser to our URL
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\nStopped.")
    finally:
        srv.server_close()
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
