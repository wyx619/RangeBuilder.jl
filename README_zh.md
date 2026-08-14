# RangeBuilder.jl

[English README](README.md)

`RangeBuilder.jl` 是对 R `alphahull` 计算接口的 JuliaGeometry 实现。几何核心不依赖
R，基于 `DelaunayTriangulation.jl`；`RCall.jl` 只用于回归与参考对照测试，不是运行时依赖。

本包还为 `DelaunayTriangulation.jl 1.6.6` 提供运行时兼容快速路径：以 O(1) 的哨兵查询
代替反复扫描整个图来判断 ghost vertex 是否存在，无须修改上游包。

## 性能

性能是本项目的首要目标。优化集中在 alpha-hull 范围构建最昂贵的两个阶段：
Delaunay/Voronoi 构建，以及动态 alpha 搜索中对候选值的反复评估。

- **一次搜索只构建一次三角剖分。** `getDynamicAlphaHull()` 只创建一次
  Delaunay/Voronoi 结构，随后所有候选 alpha 都复用该结构。搜索 *k* 个 alpha 时，
  因而只需一次三角剖分，而不是 *k* 次。
- **受保护的 Delaunay 快速路径。** 对 `DelaunayTriangulation.jl 1.6.6`，兼容层将
  重复的全图 ghost-vertex 扫描替换为 O(1) 哨兵查询，避免大点集出现严重的尺度退化。
- **alpha-hull 装配优化。** 圆弧与圆的预筛、补集索引缓存、基于端点索引的圆弧排序，
  消除了大量无效工作；点覆盖判断采用 prepared GEOS 几何，范围缓冲复用坐标转换。
- **稀疏圆弧候选。** 使用圆包围盒 sweep 构建稀疏邻接表，替代 `a × a` 的密集布尔矩阵；
  原有圆弧裁剪状态机和候选顺序保持不变。
- **严格等价的网格后处理短路。** 对大范围使用 prepared GEOS 的 `covers` 判断直接接受
  完全被范围覆盖的网格；边界网格仍使用原始 GEOS `intersects`，因此网格结果不变。

在开发机上，10,000 点、4 个 alpha 候选的热态动态搜索优化后约为 **1.5–1.6 秒**，
此前实现约为 **46.7 秒**。具体时间会随数据和硬件变化，但结构性改进是确定的：
不会再为每个 alpha 重复构建昂贵的三角剖分。

### 与原始 R 包的对照

下表是同一开发机、确定性输入、热态计时中位数的结果。它直接对比已安装的原始 R
`alphahull` 与 `rangeBuilder`，并非与另一个重实现对比。

| 工作负载 | Julia RangeBuilder.jl | 原始 R 包 | 相对结果 |
| --- | ---: | ---: | ---: |
| `ahull` + `areaahull`，500 点 | 0.0071 s | `alphahull`: 0.0600 s | 快 8.5 倍 |
| `ahull` + `areaahull`，2,000 点 | 0.0318 s | `alphahull`: 0.2200 s | 快 6.9 倍 |
| `ahull` + `areaahull`，10,000 点 | 0.2225 s | `alphahull`: 1.1300 s | 快 5.1 倍 |
| `ahull` + `areaahull`，100,000 点 | 4.9479 s | `alphahull`: 39.0300 s | 快 7.9 倍 |
| 动态搜索，500 点、4 个 alpha 候选 | 0.0453 s | `rangeBuilder`: 4.17 s | 快 92 倍 |
| 动态搜索，1,000 点、10 个 alpha 候选 | 0.4474 s | `rangeBuilder`: 30.75 s | 快 68.7 倍 |
| 动态搜索，10,000 点、10 个可用 alpha 候选 | 1.5555 s | `rangeBuilder`: 1102.1 s | 快 708 倍* |

动态搜索均使用 `fraction=0.95`、`partCount=3`、`buff=0`，并关闭海岸裁切，
以隔离范围构建本身。500 点行搜索 `0.005:0.005:0.020`；1,000 点行搜索
`0.005:0.005:0.050`。两端在这些固定点集上均得到相同的回退标签 `alphaMCH`。
大幅差异的原因在于原始 R 流程会在 alpha 改变时反复调用 `alphahull::ahull()`，
而 RangeBuilder.jl 始终复用同一个 Delaunay/Voronoi 结构。

在更大的动态搜索目标（10,000 点、10 个可用候选）中，RangeBuilder.jl 选择
`alpha0.02`，热态中位数为 **1.5555 秒**。原始 R `rangeBuilder` 的一次运行耗时
**1102.1 秒**（约 18.4 分钟），且同样选择 `alpha0.02`。因此表中的 708 倍是
单次 R 对照的指示性结果（*），不是 R 的中位数；完整 7 次 R 中位数耗时过长，
故未声称或报告。

`test/` 中的基准脚本覆盖缓存后的 alpha-hull 组件与端到端构建，并使用固定随机种子，
适合比较同一台机器上的不同修订版本。

### Rosales 全量工作流

Julia GBIF 工作流使用 8 个 Julia 线程处理 1,063,070 条清洗记录和 9,263 个物种，
成功生成 9,263 个物种范围（7,519 个 alpha-hull 任务、1,744 个稀疏物种直接任务），
并输出 763,715 行物种-网格记录。端到端耗时为 **211.46 秒**，其中范围构建耗时
**191.54 秒**。这是 Julia 工作流的实测，不是与 R 的同条件对照。

工作流脚本有意被 Git 忽略。根据本机数据路径调整后，可在 Julia 交互环境中运行：

```powershell
@'
include(raw"C:\path\to\alphahull\R\build_species_alpha_hull_1deg_ranges_gbif.jl")
using .RangeBuilderWorkflow
run_family_ranges(
    grid_path=raw"C:\path\to\alphahull\R\1d\1.shp",
    occurrence_file=raw"E:\Rosales\native_records.csv.gz",
    output_directory=raw"C:\path\to\alphahull\R\output\1deg_ranges",
    family="Rosales", workers=8, overwrite=true,
    buffer_m=10_000, fraction=0.95, part_count=3,
    initial_alpha=2, alpha_increment=1, alpha_cap=400,
    clip_to_coast=:terrestrial,
)
'@ | julia -t 8 --project=R/julia -
```

运行会在指定目录写出 `all_species_range_summary.csv`、物种网格分布文件和
`workflow_timings.csv`。

## 使用

```julia
using RangeBuilder

points = [0.0 0.0; 2.0 0.0; 2.0 1.0; 0.0 2.0; 0.6 0.7]
hull = ahull(points; alpha=0.7)
inside = inahull(hull, [0.5 0.5; 3.0 3.0])

range = getDynamicAlphaHull(points; fraction=0.95, clipToCoast=:no)
filtered = filterByProximity(points, 20.0)

ranges = (species_a=ah2polygon(ahull(points; alpha=0.7)),
          species_b=ah2polygon(ahull(points .+ 0.5; alpha=0.7)))
range_stack = rasterStackFromPolyList(ranges; resolution=0.2)
richness = speciesRichness(range_stack)

batch = buildRanges(Dict("species_a" => points, "species_b" => points[1:2, :]);
                    clipToCoast=:no)
```

`getDynamicAlphaHull()` 返回 `(hull, alpha)`：`hull` 是实现 GeoInterface 的多边形几何，
`alpha` 是最终选定的 alpha；若回退到凸包，则为 `"alphaMCH"`。`buff` 单位为米；
海岸裁切使用内置的 Natural Earth 1:50m 陆地数据，因此可离线运行。

`rasterStackFromPolyList()` 接收 GeoInterface 多边形的具名元组或字典，返回
`Rasters.RasterStack`。每个图层中，栅格中心落在多边形内的像元为 `1`，其他为
`missing`。设置 `retainSmallRanges=false` 可丢弃没有覆盖任何像元中心的小范围；
也可通过 `extent=[xmin, xmax, ymin, ymax]` 指定范围。`speciesRichness()` 将图层叠加为
每个像元的物种数；默认 `zeroToMissing=true`，无物种像元为 `missing`，与 R 示例流程一致。

`filterByLand()` 对位于经 2 km 缓冲的 Natural Earth 陆地多边形内的出现点返回 `true`，
海洋点返回 `false`，缺失坐标返回 `missing`。内置 Natural Earth v5.1.2 的 1:50m 陆地数据，
因此海岸裁切与陆地过滤均不需要网络。

`buildRanges()` 是保守的多物种入口：成功建模的多边形存入 `ranges`；少于 3 个唯一点或
几何退化的物种记入 `excluded`，不会被自动转换为推断的缓冲范围。

## 外部绘图

RangeBuilder.jl 有意不携带绘图依赖，也不导出绘图 API。这样批处理与服务器环境不会被
强制引入图形后端。它返回的多边形实现 GeoInterface 协议，丰富度结果是
`Rasters.Raster`；两者均可由调用项目自行安装的绘图库绘制。

例如，在需要地图的外部项目中安装 Makie 后端与 GeoMakie：

```julia
import Pkg
Pkg.add(["CairoMakie", "GeoMakie"])
```

随后可直接绘制计算得到的范围和原始出现点，而无须把这些包加入 RangeBuilder.jl：

```julia
using CairoMakie, GeoMakie

range = getDynamicAlphaHull(points; buff=1_000, clipToCoast=:terrestrial)
fig = Figure()
ax = GeoAxis(fig[1, 1]; dest="+proj=eqearth")
poly!(ax, range.hull; color=(:steelblue, 0.35), strokecolor=:steelblue)
scatter!(ax, points[:, 1], points[:, 2]; color=:black, markersize=7)
fig
```

如需交互式桌面或浏览器显示，可改用 `GLMakie` 或 `WGLMakie`。若只需简单查看丰富度栅格，
可以将值物化为任意绘图库可用的数组，例如用 Makie 执行 `heatmap(Array(richness))`。

## 兼容性

下列计算函数保持 R 的矩阵列模式与从 1 开始的点编号：`delvor`、`ashape`、
`complement`、`ahull`、`dw`、`inter`、`lengthahull`、`areaahulleval`、`areaahull`、
`inahull`、`rotation`、`anglesArc`、`koch` 与 `rkoch`。

`arc()` 返回采样坐标，不会向图形设备绘图。`dw_track()` 与 `ahull_track()` 同样返回
采样圆弧坐标矩阵的向量，而不是 R `ggplot2` 图层；二者都接受 `rng` 关键字以实现可重复抽样。

R 与 Julia 的 Delaunay 内部遍历顺序不同。对非退化输入，得到的几何结果等价；但
Delaunay 矩阵行顺序与插入的 alpha-hull 端点编号可能不同。

## 验证与复现

```powershell
julia --project=. -e "using Pkg; Pkg.test()"
julia --project=test test/rcall_reference.jl
julia --project=. test/benchmark_core.jl
julia --project=. test/benchmark_end_to_end.jl
Rscript test/benchmark_end_to_end.R
julia --project=test test/benchmark_dynamic_reference.jl 500 1000
```

`rcall_reference.jl` 对照已安装的 R `alphahull` 包，而不是仓库中的 R 源码。
RCall 被隔离在 `test/Project.toml`，不属于 RangeBuilder.jl 的运行时依赖。

`benchmark_core.jl` 会热身 Julia，并在固定随机种子下报告 `delvor`、`ashape`、
`complement`、`ahull` 与 `dw` 的中位时间和分配量。该基准用于同机版本间比较，
不应用于跨机器的绝对性能结论。

`benchmark_end_to_end.jl` 与 `benchmark_end_to_end.R` 使用相同的确定性点集，
比较 `ahull` 加面积计算的端到端耗时。

`benchmark_dynamic_reference.jl` 通过 RCall 对照 `getDynamicAlphaHull()` 与已安装的
原始 R `rangeBuilder`。RCall 仅用于测试环境，并非包的运行时依赖。无参数时该脚本运行
500 与 1,000 点参考情形。
