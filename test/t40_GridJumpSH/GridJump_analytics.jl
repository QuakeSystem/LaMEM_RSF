# Spurious SH-pulse reflection at a fine/coarse nseg_z jump.
# Homogeneous elastic medium; only the grid spacing jumps.

const GJ_ρ     = 3000.0
const GJ_μ     = 3.0e10
const GJ_Vs    = sqrt(GJ_μ / GJ_ρ)
const GJ_V0    = 1.0e-4
const GJ_z0    = 1000.0
const GJ_sigma = 150.0
const GJ_iface = 2000.0
const GJ_ztop  = 8000.0
const GJ_dz_f  = 50.0
const GJ_Lf    = 2000.0
const GJ_gauge = 1500.0

# Time windows at the gauge (seconds), sized for Vs and the geometry above.
# Incident: up-going half of the pulse (~0.16 s at the gauge).
# Reflected: return from the nseg jump, before the free-slip bottom bounce
# (~0.79 s) reaches the gauge.
const GJ_t_inc_lo = 0.05
const GJ_t_inc_hi = 0.28
const GJ_t_ref_lo = 0.40
const GJ_t_ref_hi = 0.62

function _replace_assign_rest(line::AbstractString, key::AbstractString, value::AbstractString)
    m = match(Regex("^(\\s*$(key)\\s*=\\s*).*\$"), line)
    m === nothing && return line
    return m.captures[1] * value
end

"""
    jump_grid(r) -> (; nseg_z, nel_z, coord_z)

Fine belt always `dz = 50` m on `0 … 2000`. Coarse belt uses `dz = r*50` on
`2000 … 8000`. `r = 1` is a uniform fine grid over the full height.
"""
function jump_grid(r::Int)
    if r == 1
        nel = Int(round(GJ_ztop / GJ_dz_f))
        return (nseg_z = 1, nel_z = string(nel), coord_z = "0.0 $(GJ_ztop)")
    end
    dz_c = r * GJ_dz_f
    nel_f = Int(round(GJ_Lf / GJ_dz_f))
    nel_c = Int(round((GJ_ztop - GJ_iface) / dz_c))
    return (
        nseg_z = 2,
        nel_z = "$nel_f $nel_c",
        coord_z = "0.0 $(GJ_iface) $(GJ_ztop)",
    )
end

function write_gridjump_dat(src::String, dst::String; jump::Int)
    g = jump_grid(jump)
    lines = String[]
    for line in eachline(src)
        line = _replace_assign_rest(line, "nseg_z", string(g.nseg_z))
        line = _replace_assign_rest(line, "nel_z", g.nel_z)
        line = _replace_assign_rest(line, "coord_z", g.coord_z)
        push!(lines, line)
    end
    write(dst, join(lines, "\n") * "\n")
    return g
end

function _vx_profile(dir::String, step; outfile::String="GridJump")
    fields = ("velocity [m/s]",)
    data, t = read_LaMEM_timestep(outfile, step, dir; fields=fields)
    tt = t isa AbstractArray ? Float64(t[1]) : Float64(t)
    z  = Float64.(data.z.val[1, 1, :])
    vx = Float64.(data.fields.velocity[1][1, 1, :])
    return z, vx, tt
end

"""
    R = gridjump_reflection(dir; outfile="GridJump")

Peak |vx| in the reflected-time window at z≈1500 m, divided by the peak |vx|
in the incident-time window at the same gauge.
"""
function gridjump_reflection(dir::String; outfile::String="GridJump")
    steps, _, _ = read_LaMEM_simulation(outfile, dir)
    isempty(steps) && error("No timesteps in $dir for $outfile")

    # gauge index from the first frame (z is fixed)
    z0, _, _ = _vx_profile(dir, steps[1]; outfile=outfile)
    ig = argmin(abs.(z0 .- GJ_gauge))

    inc_peak = 0.0
    ref_peak = 0.0
    for step in steps
        _, vx, tt = _vx_profile(dir, step; outfile=outfile)
        a = abs(vx[ig])
        if GJ_t_inc_lo <= tt <= GJ_t_inc_hi
            inc_peak = max(inc_peak, a)
        elseif GJ_t_ref_lo <= tt <= GJ_t_ref_hi
            ref_peak = max(ref_peak, a)
        end
    end
    inc_peak == 0.0 && error("No incident signal at gauge in $dir")
    return (
        jump_R = ref_peak / inc_peak,
        incident = inc_peak,
        reflected = ref_peak,
        Vs = GJ_Vs,
    )
end

"""
    res = run_gridjump_case(dir; jump)

Rewrite the input for the given cell-size jump, run LaMEM, return reflection stats.
"""
function run_gridjump_case(dir::String; jump::Int)
    cur = pwd()
    src = joinpath(cur, dir, "GridJump.dat")
    dst_name = "GridJump_run.dat"
    dst = joinpath(cur, dir, dst_name)
    write_gridjump_dat(src, dst; jump=jump)

    cd(dir)
    ok = run_lamem_local_test(dst_name, 1; outfile="GridJump.out", bin_dir="../../bin", mpiexec=mpiexec)
    cd(cur)
    ok || error("LaMEM failed for jump=$jump. See $(joinpath(dir, "GridJump.out"))")
    r = gridjump_reflection(dir; outfile="GridJump")
    return (; jump, R = r.jump_R, incident = r.incident, reflected = r.reflected)
end

"""
    cases = gridjump_sweep(dir)

Reflection coefficients for jumps 1 (uniform), 2, 3, 10 and 30.
"""
function gridjump_sweep(dir::String)
    R1  = run_gridjump_case(dir; jump=1)
    R2  = run_gridjump_case(dir; jump=2)
    R3  = run_gridjump_case(dir; jump=3)
    R10 = run_gridjump_case(dir; jump=10)
    R30 = run_gridjump_case(dir; jump=30)
    return (; R1, R2, R3, R10, R30)
end
