MESH GENERATION
===============

Objective: convert SWC nephron traces into projector-backed UGX meshes, refine
existing meshes, and inspect mesh quality.

Active entry points are directly in this folder. Supporting Python utilities
are colocated because the Lua pipelines call them by path. Earlier variants are
under archive/, and reusable repair utilities are under tools/.

See COMMANDS.txt for validated command examples and parameter explanations.
