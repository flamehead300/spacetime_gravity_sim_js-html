# spacetime_gravity_sim_js-html
Inspired by a gravity simulator I used to spend hours on, this is a browser-based N-body simulator built around the Paczyński–Wiita pseudo-Newtonian potential with optional 1PN/2PN corrections, gravitational wave energy loss, tidal disruption, and fragmentation.

# Stellar Dynamics — Gravity Engine

A real-time N-body gravity simulator with pseudo-Newtonian and post-Newtonian physics.
Two standalone files, no dependencies, no install required.

# spacetime_gravity_sim_js-html

A browser-based real-time N-body simulator built to feel like the old gravity sandboxes people used to get stuck in for hours, but with more structure and more knobs.

This version uses:
- Newtonian gravity
- Paczyński–Wiita pseudo-Newtonian gravity
- optional 1PN / boosted 2PN-style corrections
- gravitational-wave damping for heavy systems
- tidal disruption
- fragmentation
- merger ejecta
- dark matter halo background potential

Two standalone entry points. No dependencies. No install.

---

# Stellar Dynamics — Gravity Engine

Real-time N-body simulation with arcade-fast interaction and relativistic-style approximations.

## What this is

This project is meant to be:
- immediate to run
- easy to experiment with
- visually readable
- more interesting than a plain Newtonian gravity toy

It is **not** a precision GR solver. The relativistic options are approximations designed to produce useful and educational behavior in a real-time interactive sim.

---

## Files

| File | Purpose |
|---|---|
| `gravity_js.html` | Main version. Open it directly in a browser. Physics and rendering both run in JavaScript. |
| `gravity_web.py` | Python-hosted version. Starts a local server and serves the simulation in the browser. |

---

## Run

## `gravity_js.html`

Fastest way to use it.

1. Open `gravity_js.html` in Chrome, Edge, or Firefox.
2. The sim starts immediately with a binary system already spawned.

No terminal. No server. No Python.

---

## `gravity_web.py`

1. Open a terminal in the project folder.
2. Run:

```bash
python gravity_web.py

Core interaction
Canvas controls
Action	Result
Click and release	Spawn a body with zero initial velocity
Click, drag, release	Spawn a body with launch velocity opposite the drag direction
Scroll / trackpad / pinch	Zoom toward the cursor

Longer drag = faster launch.

The camera follows the first body in the internal body list, which usually means the initial seeded body unless bodies merge, fragment, or get removed.

UI panels
BODY PARAMS
SPAWN MASS

Sets the mass of the next body you place.

The panel also shows:

current visual radius

Schwarzschild radius estimate r_s = 2Gm / c²

PHYSICS ENGINE
Control	Meaning
GRAVITY [G]	Global gravity strength
SOFTENING [e]	Softening term added to short-range interactions to suppress force spikes
TIME STEP [dt]	Base simulation timestep
DRAG COEFF	Velocity damping multiplier
MAX VELOCITY	Hard speed cap
Notes

Higher G pulls everything together faster

Higher dt makes the sim advance faster but can destabilize close encounters

DRAG COEFF < 1 bleeds momentum over time

MAX VELOCITY helps prevent runaway explosions

RELATIVITY

These are approximations, not full GR.

Control	Meaning
ENABLE RELATIVITY	Switches from pure Newtonian gravity to relativistic-style approximations
1PN / EIH CORRECTIONS	Adds first post-Newtonian style corrections
2PN HIGHER-ORDER TERMS	Adds an extra boosted higher-order term on top of 1PN
GW DAMPING (HEAVY)	Adds strong gravitational-wave-like orbital energy loss for sufficiently massive systems
SHOW EVENT HORIZONS	Draws dashed horizon circles
LIGHTSPEED [c]	Effective simulation light speed
Mode behavior
Relativity off

Standard softened Newtonian gravity.

Relativity on, 1PN off

Uses a Paczyński–Wiita-style pseudo-Newtonian potential. This makes gravity steepen near the Schwarzschild radius and gives black-hole-like behavior without solving Einstein field equations.

1PN on

Adds post-Newtonian style correction terms. This is where you start seeing more interesting close-orbit behavior like precession.

2PN on

Enables an additional higher-order boost on top of 1PN. This is an approximation layer for stronger compact-object behavior, not a rigorous full 2PN expansion.

Event horizons

When enabled, the sim draws a dashed purple circle using:

r_s = 2Gm / c²

These are visualized horizon estimates for the current sim parameters.

COLLISION
Control	Meaning
FRAGMENTATION	High-speed impacts can shatter bodies into fragments instead of merging
TIDAL DISRUPTION	Small bodies crossing within the Roche-style disruption limit of a heavier body break apart
MERGER EJECTA	Energetic mergers can throw off a small debris body
RESTITUTION [e]	Collision elasticity
FRAG SPEED	Velocity threshold used to trigger fragmentation
Collision modes
RESTITUTION = 0

Bodies merge on impact unless fragmentation logic overrides.

RESTITUTION > 0

Bodies bounce with partial or full elasticity depending on the value.

Fragmentation

If enabled and the collision speed passes the threshold, the colliding pair can be replaced with fragments.

Tidal disruption

If enabled, a lighter body can be torn into pieces before direct collision if it crosses the disruption radius of a heavier one.

ADVANCED
Control	Meaning
ADAPTIVE DT	Shrinks timestep in high-acceleration situations
LEAPFROG INTEGRATOR	Uses leapfrog-style integration instead of the simpler update path
REALISTIC RADII	Changes how mass maps to radius
DARK MATTER HALO	Adds a background Plummer-like halo potential centered on the origin
Notes
Adaptive dt

Useful for close encounters and dense systems where fixed timesteps become unstable.

Leapfrog integrator

Usually gives better orbital behavior than the simpler direct update path, especially when drag is low and you care about longer-term motion.

Realistic radii

Switches away from a simple visual cube-root mass scaling and uses a more stylized regime-based mapping.

Dark matter halo

Adds a smooth central background potential that can help keep wider systems loosely bound.

REALITY CONTROLS
Button	Meaning
START / HOLD ENGINES	Toggle simulation run state
SINGLE STEP	Advance one simulation step
CLEAR	Remove all bodies
SPAWN 100	Scatter 100 random bodies into the visible region
Spawn 100 behavior

Bulk spawn temporarily pauses the sim, adds bodies in a chunk, updates the display, then resumes if the sim was already running.

TELEMETRY

The telemetry panel shows:

body count

render FPS

physics time per step

current performance mode

Modes
Mode	Meaning
NORMAL	Full stepping is active
THROTTLED	Body count is high enough that per-frame stepping is reduced or skipped
BODY CAP	Maximum body count reached
Performance model

This sim is still fundamentally O(n²) for force evaluation.

That means:

double the bodies -> roughly four times the pair work

collision checking is also O(n²)

Current throttle behavior
Bodies	Behavior
<= 80	4 physics steps per frame
81–150	2 physics steps per frame
151–300	1 physics step per frame
> 300	physics stepping is throttled hard
>= 1200	hard body cap
Collision cutoff

Collision resolution is automatically disabled above:

250 bodies

This is intentional. It keeps the browser responsive instead of letting overlap resolution dominate the frame.

Physics notes
Force model

At the base level, forces are pairwise and softened:

gravity acts between every body pair

softening prevents singular short-range spikes

optional halo term adds a global background potential

Relativistic options

These are there to create richer orbital behavior in a fast interactive sim.

They should be read as:

pseudo-Newtonian / post-Newtonian inspired

educational / demonstrative

not a strict research-grade relativistic integrator

Collision logic

Bodies can:

merge

bounce

fragment

undergo tidal disruption

emit ejecta in energetic mergers

That combination makes dense scenes much more interesting than simple point-mass gravity alone.

Troubleshooting
Browser becomes slow or stalls

Most likely causes:

too many bodies

too many close interactions

collision-heavy scene before the collision cutoff kicks in

Fix

hit CLEAR

restart with fewer bodies

reduce mass or reduce fragmentation-heavy scenarios

Nothing is moving

Check these first:

sim may be paused

TIME STEP [dt] may be very low

DRAG COEFF may be draining motion quickly

body count may be high enough that physics is throttled

Bodies instantly merge

Likely causes:

RESTITUTION = 0

bodies were spawned overlapping

softening and initial placement created rapid collapse

Fix

increase separation

raise restitution

lower spawn mass

increase softening slightly

Event horizons are not visible

Requirements:

ENABLE RELATIVITY must be on

SHOW EVENT HORIZONS must be on

Also note:

low-mass bodies may have very small horizons

you may need to zoom in or raise mass

Orbits blow up immediately

Most likely:

G is too high for the initial setup

dt is too large

c is too low and relativistic terms are exaggerated

launch speed was too high

Fix

Reset toward:

moderate G

lower dt

higher c

lower spawn velocities

gravity_web.py will not start
Port already in use

Another process is already bound to 127.0.0.1:8765.

Fix

Kill the old process or change the port in the Python file.

Blank page after browser opens

The Python process may have crashed after opening the browser.

Fix

Check the terminal output for the traceback, then refresh once the server is running again.

Wrong Python version

Use:

python --version

You need Python 3.9+.

Suggested experiments
Stable binary

Spawn two heavy bodies at moderate separation with tangential velocity.

Inspiral

Turn on:

relativity

GW damping

Then use two very massive bodies and lower c enough to make the effect visible.

Tidal shredding

Use one very heavy body and one much lighter one. Place the lighter one on a close pass.

Fragmentation test

Turn on fragmentation, lower the threshold, and throw a fast body directly into another.

Halo confinement

Enable dark matter halo and scatter many bodies with mild velocities to create a looser bound cloud.

Summary

This project is a fast, self-contained gravity sandbox with:

real-time interaction

zero dependencies

relativistic-style approximations

collision / fragmentation systems

performance throttling for browser survival

Open the HTML file and start throwing bodies around.


If you want, I can also do a second pass that makes it sound even mor
