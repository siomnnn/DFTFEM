using Plots

#
#point1 = [-0.4233581327227974, 0.5, 0.1499999999990289]
#point2 = [-0.5, 0.4116880127287288, 0.05233514768171421]
#point3 = [-0.5, 0.413397459352698, 0.1499999996330086]
#point4 = [-0.336395025265869, 0.3336207719605636, 0.1348284560558868]
#pointx = [-0.4197584844840181, 0.46480253970766194, 0.1191254038545938]
#
#plot(dpi=400, aspect_ratio=:equal)
#scatter!([-0.4233581327227974], [0.5], [0.1499999999990289], color=:green)
#scatter!([-0.5], [0.4116880127287288], [0.05233514768171421], color=:green)
#scatter!([-0.5], [0.413397459352698], [0.1499999996330086], color=:green)
#scatter!([-0.336395025265869], [0.3336207719605636], [0.1348284560558868], color=:green)
#scatter!([-0.4197584844840181], [0.46480253970766194], [0.1191254038545938], color=:red)
#
#println([-0.4233581327227974 - 0.5 - 0.5 - 0.336395025265869, 0.5 + 0.4116880127287288 + 0.413397459352698 + 0.3336207719605636, 0.1499999999990289 + 0.05233514768171421 + 0.1499999996330086 + 0.1348284560558868]/4)
#
#println(norm(point1-pointx))
#println(norm(point2-pointx))
#println(norm(point3-pointx))
#println(norm(point4-pointx))
#
#
#xlims!(-0.5, 0.5)
#ylims!(-0.5, 0.5)
#zlims!(-0.5, 0.5)

#point0 = [0.1381966011250105, 0.1381966011250105, 0.1381966011250105]
#point1 = [0.5854101966249685, 0.1381966011250105, 0.1381966011250105]
#point2 = [0.1381966011250105, 0.5854101966249685, 0.1381966011250105]
#point3 = [0.1381966011250105, 0.1381966011250105, 0.5854101966249685]
#
#point0 = VectorValue{3, Float64}(3.7979370893331734, 0.008273567508750156, -0.5925389225905165)
#point1 = VectorValue{3, Float64}(3.7979370893307447, -0.45070633319845305, -0.1335590218857422)
#point2 = VectorValue{3, Float64}(3.7979370893332924, 0.46725346821328756, -0.13355902188586088)
#point3 = VectorValue{3, Float64}(2.803480637425081, 0.03504739438616499, -0.1067851950114933)
#
#plot()
#scatter!([point0[1]], [point0[2]], [point0[3]], color=:red, label="point0")
#scatter!([point1[1]], [point1[2]], [point1[3]], color=:green, label="point1")
#scatter!([point2[1]], [point2[2]], [point2[3]], color=:blue, label="point2")
#scatter!([point3[1]], [point3[2]], [point3[3]], color=:yellow, label="point3")
#
#cell_map = get_cell_map(Ω)[1]
#
#corner0 = cell_map(VectorValue(0, 0, 0))
#corner1 = cell_map(VectorValue(1, 0, 0))
#corner2 = cell_map(VectorValue(0, 1, 0))
#corner3 = cell_map(VectorValue(0, 0, 1))
#
#scatter!([corner0[1]], [corner0[2]], [corner0[3]], color=:black)
#scatter!([corner1[1]], [corner1[2]], [corner1[3]], color=:black)
#scatter!([corner2[1]], [corner2[2]], [corner2[3]], color=:black)
#scatter!([corner3[1]], [corner3[2]], [corner3[3]], color=:black)

energies = [43.80033057191101, 36.46090934678961, 32.50361011140525, 31.252795225124224, 30.78166277853551, 30.410500755548846, 30.179750692134576]
plot([4, 6, 8, 10, 12, 14, 16], abs.(energies .- 30.09222653256), xlabel="N", ylabel="E_{FEM} - E_{DFTK} (Hartree)", title="error vs N", legend=false, xscale=:log10, yscale=:log10)
savefig("julia/energy_vs_N.png")