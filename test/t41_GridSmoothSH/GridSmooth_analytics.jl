# Overall fine→coarse factor 50 for an elastic SH pulse.
# Compare an abrupt jump, a biased smooth segment, jump-5 segments, and
# 10 geometric segments (max local ratio ~1.54).

const GS_ρ     = 3000.0
const GS_μ     = 3.0e10
const GS_Vs    = sqrt(GS_μ / GS_ρ)
const GS_dz_f  = 50.0
const GS_Lf    = 2000.0
const GS_gauge = 1500.0
const GS_ratio = 50.0

# Incident / reflected windows at the fine-side gauge (same geometry as t40).
const GS_t_inc_lo = 0.05
const GS_t_inc_hi = 0.28
# Stop before the free-slip bottom bounce (~0.79 s) reaches the gauge.
const GS_t_ref_lo = 0.40
const GS_t_ref_hi = 0.58

function _replace_assign_rest(line::AbstractString, key::AbstractString, value::AbstractString)
    m = match(Regex("^(\\s*$(key)\\s*=\\s*).*\$"), line)
    m === nothing && return line
    return m.captures[1] * value
end

function _replace_layer_top(line::AbstractString, ztop::Float64)
    m = match(r"^(\s*top\s*=\s*)\S+", line)
    m === nothing && return line
    return m.captures[1] * string(ztop)
end

"""
    grid_case(name) -> NamedTuple

Return nseg_z, nel_z, coord_z, bias_z and ztop for one expansion strategy.
"""
function grid_case(name::AbstractString)
    dz_f = GS_dz_f
    Lf   = GS_Lf
    nel_f = Int(round(Lf / dz_f))
    dz_c = GS_ratio * dz_f   # 2500 m

    if name == "abrupt"
        # one interface: 50 → 2500
        nel_c = 20
        ztop = Lf + nel_c * dz_c
        return (
            name = name,
            nseg_z = 2,
            nel_z = "$nel_f $nel_c",
            coord_z = "0.0 $Lf $ztop",
            bias_z = "1.0 1.0",
            ztop = ztop,
        )
    elseif name == "bias"
        # begSz=50, endSz=2500, bias=50 → avg=1275; M=40 → max neighbour ~2.26
        M = 40
        avg = dz_f * (1.0 + GS_ratio) / 2.0
        Ltr = avg * M
        ztop = Lf + Ltr
        return (
            name = name,
            nseg_z = 2,
            nel_z = "$nel_f $M",
            coord_z = "0.0 $Lf $ztop",
            bias_z = "1.0 $(GS_ratio)",
            ztop = ztop,
        )
    elseif name == "segs5"
        # uniform segments: ×5, ×5, ×2  (50 → 250 → 1250 → 2500)
        dzz = (250.0, 1250.0, 2500.0)
        nels = [nel_f]
        coords = [0.0, Lf]
        z = Lf
        for dz in dzz
            n = 8
            z += n * dz
            push!(nels, n)
            push!(coords, z)
        end
        nb = length(nels)
        return (
            name = name,
            nseg_z = nb,
            nel_z = join(nels, " "),
            coord_z = join(string.(coords), " "),
            bias_z = join(fill("1.0", nb), " "),
            ztop = z,
        )
    elseif name == "segs10"
        # 10 segments: fine pad + 9 geometric steps, ratio 50^(1/9) ≈ 1.54
        nels = [nel_f]
        coords = [0.0, Lf]
        z = Lf
        for k in 1:9
            dz = dz_f * GS_ratio^(k / 9)
            n = 4
            z += n * dz
            push!(nels, n)
            push!(coords, z)
        end
        nb = length(nels)
        @assert nb <= 10
        return (
            name = name,
            nseg_z = nb,
            nel_z = join(nels, " "),
            coord_z = join(string.(coords), " "),
            bias_z = join(fill("1.0", nb), " "),
            ztop = z,
        )
    else
        error("Unknown grid case $name")
    end
end

function write_gridsmooth_dat(src::String, dst::String; name::AbstractString)
    g = grid_case(name)
    lines = String[]
    for line in eachline(src)
        line = _replace_assign_rest(line, "nseg_z", string(g.nseg_z))
        line = _replace_assign_rest(line, "nel_z", g.nel_z)
        line = _replace_assign_rest(line, "coord_z", g.coord_z)
        line = _replace_assign_rest(line, "bias_z", g.bias_z)
        line = _replace_layer_top(line, g.ztop)
        push!(lines, line)
    end
    write(dst, join(lines, "\n") * "\n")
    return g
end

function _vx_profile(dir::String, step; outfile::String="GridSmooth")
    fields = ("velocity [m/s]",)
    data, t = read_LaMEM_timestep(outfile, step, dir; fields=fields)
    tt = t isa AbstractArray ? Float64(t[1]) : Float64(t)
    z  = Float64.(data.z.val[1, 1, :])
    vx = Float64.(data.fields.velocity[1][1, 1, :])
    return z, vx, tt
end

function gridsmooth_reflection(dir::String; outfile::String="GridSmooth")
    steps, _, _ = read_LaMEM_simulation(outfile, dir)
    isempty(steps) && error("No timesteps in $dir for $outfile")

    z0, _, _ = _vx_profile(dir, steps[1]; outfile=outfile)
    ig = argmin(abs.(z0 .- GS_gauge))

    inc_peak = 0.0
    ref_peak = 0.0
    for step in steps
        _, vx, tt = _vx_profile(dir, step; outfile=outfile)
        a = abs(vx[ig])
        if GS_t_inc_lo <= tt <= GS_t_inc_hi
            inc_peak = max(inc_peak, a)
        elseif GS_t_ref_lo <= tt <= GS_t_ref_hi
            ref_peak = max(ref_peak, a)
        end
    end
    inc_peak == 0.0 && error("No incident signal at gauge in $dir")
    return (jump_R = ref_peak / inc_peak, incident = inc_peak, reflected = ref_peak)
end

function run_gridsmooth_case(dir::String; name::AbstractString)
    cur = pwd()
    src = joinpath(cur, dir, "GridSmooth.dat")
    dst_name = "GridSmooth_run.dat"
    dst = joinpath(cur, dir, dst_name)
    g = write_gridsmooth_dat(src, dst; name=name)

    cd(dir)
    ok = run_lamem_local_test(dst_name, 1; outfile="GridSmooth.out", bin_dir="../../bin", mpiexec=mpiexec)
    cd(cur)
    ok || error("LaMEM failed for case=$name. See $(joinpath(dir, "GridSmooth.out"))")
    r = gridsmooth_reflection(dir; outfile="GridSmooth")
    return (; name, R = r.jump_R, incident = r.incident, reflected = r.reflected, ztop = g.ztop)
end

"""
    cases = gridsmooth_sweep(dir)

Reflection coefficients for abrupt / bias / segs5 / segs10 expansions to ×50.
"""
function gridsmooth_sweep(dir::String)
    abrupt = run_gridsmooth_case(dir; name="abrupt")
    bias   = run_gridsmooth_case(dir; name="bias")
    segs5  = run_gridsmooth_case(dir; name="segs5")
    segs10 = run_gridsmooth_case(dir; name="segs10")
    return (; abrupt, bias, segs5, segs10)
end
