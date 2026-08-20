# RangeBuilder.jl

RangeBuilder.jl 是用 JuliaGeometry 实现的高性能 R `alphahull` 与 `rangeBuilder` 计算接口。几何核心完全独立于 R，运行时不需要 `RCall`。

## 安装

```julia
using Pkg
Pkg.add(url="https://github.com/wyx619/RangeBuilder.jl")
```

本地开发：

```julia
using Pkg
Pkg.develop(path="C:/path/to/RangeBuilder.jl")
```

## 快速开始

```julia
using RangeBuilder

points = [
    116.0 39.0
    116.5 39.0
    116.4 39.4
    116.1 39.3
]

hull = ahull(points; alpha=0.7)
polygon = ah2polygon(hull)
inside = inahull(hull, [116.2 39.2; 118.0 40.0])

range = getDynamicAlphaHull(
    points;
    fraction=0.95,
    partCount=3,
    buff=1_000,
    clipToCoast=:terrestrial,
)

range.alpha
range.hull
```

输入是 `n x 2` 矩阵，列顺序为 `(longitude, latitude)`。`buff` 的单位是米。返回值包含实现 GeoInterface 的多边形或多多边形，以及最终选择的 alpha 标签。

## 性能优化

主要优化集中在 Delaunay/Voronoi 构建和动态 alpha 评估：

- 缓存后端只构建一次 Delaunay 结构，所有候选 alpha 复用它；
- 为 `DelaunayTriangulation.jl 1.6.6` 提供 O(1) ghost-vertex 哨兵查询兼容路径；
- alpha-hull 装配使用圆弧预筛、补集索引缓存和基于端点的排序；
- 圆弧相交候选使用稀疏包围盒邻接表，而不是密集的 `a x a` 矩阵；
- prepared GEOS 谓词加速重复覆盖和网格边界检查，同时保持原有边界语义。

需要严格复刻 R `rangeBuilder` 候选调用顺序时使用 `backend=:shull`。它会重建每个候选，速度较慢，但保留候选级重试和 RNG 语义。

开发机上，10,000 个点、4 个候选 alpha 的热态动态搜索约为 1.5-1.6 秒，早期实现约为 46.7 秒。实际速度取决于硬件和输入几何。

## 动态范围 API

```julia
range = getDynamicAlphaHull(
    points;
    fraction=0.95,
    partCount=3,
    initialAlpha=2,
    alphaIncrement=1,
    alphaCap=400,
    buff=10_000,
    clipToCoast=:terrestrial,
    backend=:shull,
)
```

`clipToCoast` 可取 `:no`、`:terrestrial` 或 `:aquatic`。单个物种内部的动态搜索保持有序；应用程序可以使用 Julia 线程调度相互独立的物种。要严格复刻 R 的重试过程，请传入匹配的 `RSeed` RNG。

## S2 海岸语义

接受的范围先在 Equal Earth 中缓冲，再转换到 WGS84 执行最终地理海岸 overlay。`S2Geography.jl` 使用 Google S2 语义执行陆地或海洋裁剪，可处理大圆弧边界、反经线和多组件结果，运行时不需要 R 或 `RCall`。

S2 只用于最终海岸裁剪。高频的 alpha 有效性与出现点覆盖率检查使用更轻量的 GeometryOps 球面谓词；MCH 回退使用专用的球面凸包实现。

## 内置空间数据

包内置 Natural Earth 1:50m 陆地数据，用于离线海岸裁剪和陆地点过滤。可选的范围构建工作流所需的一度网格也以内部 JLD2 资源随包提供。资源路径由包自动解析，调用者不需要查找或传入文件路径。

## 多物种与物种丰富度

```julia
ranges = (
    species_a=ah2polygon(ahull(points; alpha=0.7)),
    species_b=ah2polygon(ahull(points .+ [0.5 0.5]; alpha=0.7)),
)

stack = rasterStackFromPolyList(ranges; resolution=0.2)
richness = speciesRichness(stack)
```

`rasterStackFromPolyList` 返回 `Rasters.RasterStack`；`speciesRichness` 将其转换为逐网格物种数。包本身不包含绘图后端。

## 外部绘图

多边形输出实现 GeoInterface，可以交给应用程序自己的绘图栈。例如，应用程序可以额外安装 `CairoMakie` 和 `GeoMakie`，再将 `range.hull` 传给 `poly!`。绘图保持在包外，避免批处理和服务器环境被迫安装图形后端。

## 兼容性

`delvor`、`ashape`、`complement`、`ahull`、`dw`、`inter`、`lengthahull`、`areaahulleval`、`areaahull`、`inahull`、`rotation`、`anglesArc`、`koch` 和 `rkoch` 保持 R 风格的矩阵结构与从 1 开始的点编号。

对普通输入，R 与 Julia 会产生等价几何。共圆或近退化点集可能存在多个合法 Delaunay 三角剖分，因此内部三角形、圆弧和矩阵行序不属于兼容性契约；兼容目标是最终接受的几何和范围语义。

## 开发检查

```julia
using Pkg
Pkg.test()
```


