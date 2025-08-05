using Test
using ..Yggdrasil

# Helper constructors for tests
function make_patch(id, llc, size, level)
    NBOR = NBOR_NML(0, Int[])
    return Patch_NML(id, llc .+ size ./ 2, size, level, 0.0, 0, size, [1,1,1], [1,1,1],
                     0, 0, zeros(3), 0.0, 0, "", "", 0, 0, [0,0,0], 0.0,
                     llc .+ size ./ 2, llc, zeros(3), zeros(3), zeros(3), NBOR, 0, "")
end

function make_snapshot(patches)
    io = IO_NML(0,0,0.0,false,0,"",0,false,false,false)
    snap = SNAPSHOT_NML(0,0,0.0,0, zeros(3), Int[], Int[], Int[], Int[], Int[], 0,0,0,0.0,"","", Bool[], false,0,false,0,0,0, Int[], 1, zeros(3))
    idx = IDX_NML(fill(0,27)...)
    scaling = SCALING_NML(fill(0.0,14)...,"",0,0,0.0)
    return Snapshot_metadata(io, snap, idx, scaling, length(patches), 0, patches,
                             "", false, false, 0, 0, Int[], 0, "", Particles_NML[], "", 0, 0)
end

patch1 = make_patch(1, [0.0,0.0,0.0], [1.0,1.0,1.0], 0)
patch2 = make_patch(2, [1.0,0.0,0.0], [1.0,1.0,1.0], 1)
patch3 = make_patch(3, [0.0,1.0,0.0], [1.0,1.0,1.0], 1)
meta = make_snapshot([patch1, patch2, patch3])

@testset "square and cube patch search" begin
    square = [[0.5,0.5,0.5], [1.5,0.5,0.5], [1.5,1.5,0.5], [0.5,1.5,0.5]]
    patches = find_patches_square(meta, square; all_levels=true)
    @test length(patches) == 3
    patches_lvl = find_patches_square(meta, square; level=1)
    @test Set(p.ID for p in patches_lvl) == Set([2,3])

    cube = [[1.2,0.2,0.2],[1.8,0.2,0.2],[1.8,0.8,0.2],[1.2,0.8,0.2],
            [1.2,0.2,0.8],[1.8,0.2,0.8],[1.8,0.8,0.8],[1.2,0.8,0.8]]
    cube_patches = find_patches_cube(meta, cube; all_levels=true)
    @test length(cube_patches) == 1 && cube_patches[1].ID == 2
    cube_lvl = find_patches_cube(meta, cube; level=1)
    @test length(cube_lvl) == 1 && cube_lvl[1].ID == 2
end

@testset "line patch search" begin
    startp = [-0.5,0.5,0.5]
    endp   = [1.2,0.5,0.5]
    line_patches = find_patches_line(meta, startp, endp; all_levels=true)
    @test Set(p.ID for p in line_patches) == Set([1,2])
    line_lvl = find_patches_line(meta, startp, endp; level=1)
    @test length(line_lvl) == 1 && line_lvl[1].ID == 2
end
