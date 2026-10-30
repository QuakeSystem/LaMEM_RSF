# Homogeneous elastic standing SH wave.
# vx(z,t) = V0 sin(k z) cos(ω t),  k = n π / H,  ω = sqrt(G/ρ) k
# Same numbers as SHWave.dat (SI, unit scales of 1).

const SH_ρ  = 3000.0
const SH_μ  = 3.0e10
const SH_H  = 1000.0
const SH_V0 = 1.0e-4
const SH_n  = 1
const SH_Vs = sqrt(SH_μ / SH_ρ)
const SH_k  = SH_n * π / SH_H
const SH_ω  = SH_Vs * SH_k
const SH_T  = 2π / SH_ω

"Amplitude of τ_xz."
shwave_tau_amp() = SH_μ * SH_k * SH_V0 / SH_ω

"Amplitude of ux."
shwave_disp_amp() = SH_V0 / SH_ω

function shwave_fields(z, t)
    s = sin.(SH_k .* z)
    c = cos.(SH_k .* z)
    vx = SH_V0 .* s .* cos(SH_ω * t)
    ux = shwave_disp_amp() .* s .* sin(SH_ω * t)
    tx = shwave_tau_amp() .* c .* sin(SH_ω * t)
    return vx, ux, tx
end

"Elastic plus kinetic energy per unit horizontal area, from corner samples."
function shwave_energy(vx, txz, z)
    e = 0.5 * SH_ρ .* vx .^ 2 .+ txz .^ 2 ./ (2 * SH_μ)
    dz = diff(z)
    return sum(0.5 .* (e[1:end-1] .+ e[2:end]) .* dz)
end

function shwave_rms(a)
    return sqrt(sum(abs2, a) / length(a))
end

"""
    rel = shwave_errors(dir, outfile="SHWave")

Relative RMS errors of vx, ux and τ_xz at the last output, scaled by their
analytical amplitudes, plus the relative drift of elastic-plus-kinetic energy
from the initial output.
"""
function shwave_errors(dir::String, outfile::String="SHWave")
    steps, _, _ = read_LaMEM_simulation(outfile, dir)
    fields = ("velocity [m/s]", "dev_stress [Pa]", "tot_displ [m]")
    data0, _ = read_LaMEM_timestep(outfile, 0, dir; fields=fields)
    data, t = read_LaMEM_timestep(outfile, steps[end], dir; fields=fields)
    tt = t isa AbstractArray ? Float64(t[1]) : Float64(t)

    z  = Float64.(data.z.val[1, 1, :])
    vx = Float64.(data.fields.velocity[1][1, 1, :])
    ux = Float64.(data.fields.tot_displ[1][1, 1, :])
    # dev_stress component 3 is τ_xz (ParaView order xx, xy, xz, ...)
    tx = Float64.(data.fields.dev_stress[3][1, 1, :])

    z0  = Float64.(data0.z.val[1, 1, :])
    vx0 = Float64.(data0.fields.velocity[1][1, 1, :])
    tx0 = Float64.(data0.fields.dev_stress[3][1, 1, :])

    vxa, uxa, txa = shwave_fields(z, tt)
    return (
        t = tt,
        vx = shwave_rms(vx .- vxa) / SH_V0,
        ux = shwave_rms(ux .- uxa) / shwave_disp_amp(),
        tx = shwave_rms(tx .- txa) / shwave_tau_amp(),
        energy = abs(shwave_energy(vx, tx, z) - shwave_energy(vx0, tx0, z0)) /
                 shwave_energy(vx0, tx0, z0),
    )
end

"""
    write_shwave_dat(src, dst; nel_z, nstep)

Copy the reference input and replace the grid and the fixed time step.
`nstep` steps cover `SH_T/8`.
"""
function _replace_assign(line::AbstractString, key::AbstractString, value::AbstractString)
    m = match(Regex("^(\\s*$(key)\\s*=\\s*)\\S+"), line)
    m === nothing && return line
    return m.captures[1] * value
end

"Round-trip decimal, so the rewritten input does not need Printf."
fmt_sci(x) = sprint(show, Float64(x))

function write_shwave_dat(src::String, dst::String; nel_z::Int, nstep::Int)
    dt = (SH_T / 8) / nstep
    tend = nstep * dt
    dt_s = fmt_sci(dt)
    lines = String[]
    for line in eachline(src)
        line = _replace_assign(line, "nel_z", string(nel_z))
        line = _replace_assign(line, "time_end", fmt_sci(tend))
        line = _replace_assign(line, "dt", dt_s)
        line = _replace_assign(line, "dt_min", dt_s)
        line = _replace_assign(line, "dt_max", dt_s)
        line = _replace_assign(line, "nstep_max", string(nstep + 5))
        push!(lines, line)
    end
    write(dst, join(lines, "\n") * "\n")
    return dt, tend
end

"""
    err = run_shwave_case(dir; nel_z, nstep)

Run one standing-wave case and return the error named tuple from `shwave_errors`.
"""
function run_shwave_case(dir::String; nel_z::Int, nstep::Int)
    cur = pwd()
    src = joinpath(cur, dir, "SHWave.dat")
    dst_name = "SHWave_run.dat"
    dst = joinpath(cur, dir, dst_name)
    write_shwave_dat(src, dst; nel_z=nel_z, nstep=nstep)

    cd(dir)
    ok = run_lamem_local_test(dst_name, 1; outfile="SHWave.out", bin_dir="../../bin", mpiexec=mpiexec)
    cd(cur)
    ok || error("LaMEM failed for nel_z=$nel_z, nstep=$nstep. See $(joinpath(dir, "SHWave.out"))")
    return shwave_errors(dir, "SHWave")
end

"""
    cases = shwave_convergence(dir)

Four runs over one eighth of a period:

- `coarse`: nel_z =  8, 32 steps
- `space`:  nel_z = 16, 32 steps (same dt)
- `fine`:   nel_z = 32, 32 steps
- `time`:   nel_z = 32,  8 steps (dt four times larger)
"""
function shwave_convergence(dir::String)
    coarse = run_shwave_case(dir; nel_z=8,  nstep=32)
    space  = run_shwave_case(dir; nel_z=16, nstep=32)
    fine   = run_shwave_case(dir; nel_z=32, nstep=32)
    time   = run_shwave_case(dir; nel_z=32, nstep=8)
    return (; coarse, space, fine, time)
end
