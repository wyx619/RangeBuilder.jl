# RangeBuilder.jl

RangeBuilder.jl 是高性能、与 R 兼容的 `alphahull` 与 `rangeBuilder` Julia 实现。alpha-hull 核心是对 R `interp::tri.mesh()` SHull 的纯 Julia 移植，运行时不需要 R 或 `RCall`。

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

所有 alpha-hull 调用都使用唯一、确定性的 SHull 路径。动态 alpha 搜索会在同一次
搜索内复用首次成功且未 jitter 的 SHull 三角剖分；重复点恢复、jitter 重试、拓扑失败和
MCH 回退仍保持逐候选路径，从而兼容 R `rangeBuilder` 的语义和 RNG 消耗。

范围工作流则针对其余空间热点优化：

- CEA prepared-land `covers` 会让完全在陆地内的范围跳过 S2 overlay；
- 海岸范围先查询线程本地陆地组件 `STRtree`，再进行局部 S2 overlay；不满足安全条件时自动回退全球 S2；
- worker 复用 PROJ、GEOS、S2 与网格空间索引资源；
- 已删除未使用的全球网格质心构建。

在 Rosales 基准输入（1,063,070 条记录、9,263 个物种）上，内置的
Tables.jl 工作流使用 8 个 Julia worker 成功完成 **9,263/9,263** 个物种：
`build_ranges=101.818 秒`，总耗时 `124.404 秒`。实际速度取决于硬件、输入
几何、海岸邻近程度和 worker 数。

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
)
```

`clipToCoast` 可取 `:no`、`:terrestrial` 或 `:aquatic`。单个物种内部的动态搜索保持有序；应用程序可以使用 Julia 线程调度相互独立的物种。要严格复刻 R 的重试过程，请传入匹配的 `RSeed` RNG。

## S2 海岸语义

接受的范围先在 Equal Earth 中缓冲，再转换到 WGS84 执行最终地理海岸 overlay。`S2Geography.jl` 使用 Google S2 语义执行陆地或海洋裁剪，可处理大圆弧边界、反经线和多组件结果，运行时不需要 R 或 `RCall`。

S2 只用于最终海岸裁剪。高频的 alpha 有效性与出现点覆盖率检查使用更轻量的 GeometryOps 球面谓词；MCH 回退使用专用的球面凸包实现。

## 内置空间数据

包内置 Natural Earth 1:50m 陆地数据，用于离线海岸裁剪和陆地点过滤。可选的范围构建工作流所需的一度网格也以内部 JLD2 资源随包提供。资源路径由包自动解析，调用者不需要查找或传入文件路径。

## 多物种与物种丰富度

包内还提供不依赖文件格式的内置多物种范围工作流，输入可以是任何
兼容 Tables.jl 的观测表：

```julia
context = prepare_workflow_context()
input = standardize_occurrences(
    df1;
    species=:Accepted_species,
    longitude=:mylongitude,
    latitude=:mylatitude,
    family=:Accepted_family,
    order=nothing,
)
result = build_species_ranges(input, context; workers=8)
```

`standardize_occurrences` 会先对经纬度四舍五入，再进行去重。默认保留小数点后
`2` 位；可通过 `coordinate_digits=1`（或其他非负整数）调整精度。去重键为
`(species, rounded_longitude, rounded_latitude)`，后续范围计算也使用四舍五入后的
坐标。`source_rows` 保存每个去重坐标组合首次保留记录在原表中的行号。

它使用包内置的一度网格和 Natural Earth 陆地数据，并启用 CEA
`covers` 快速路径、区域 S2 候选缓存和安全的全球 S2 回退。详见
[`src/workflow/README.md`](src/workflow/README.md)。

每个 `SpeciesRangeResult` 都会保留标准化输入中折叠后的元数据和原始行号。
工作流不负责文件读写；文件读取和结果写出应由应用层完成。

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

包使用与 R 兼容的 SHull 三角剖分路径。对普通输入，R 与 Julia 会产生等价几何；重复或近退化输入则遵循 R 风格的点恢复与 jitter 处理。兼容目标是最终接受的几何和范围语义。

## 开发检查

```julia
using Pkg
Pkg.test()
```
