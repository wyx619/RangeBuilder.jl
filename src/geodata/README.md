# Bundled geographic data

`ne_50m_land.jld2` stores the dissolved WKB geometry for Natural Earth
v5.1.2's 1:50m physical land dataset. `RangeBuilder` loads it on demand for
land filtering and coastline clipping, without a network request.

Natural Earth data is public domain. Source:
https://github.com/nvkelso/natural-earth-vector/tree/v5.1.2

Original dataset SHA-256:
`e874b27a51d146452be360cafb3cc50c86001074a67d534113e6534682f9826b`
