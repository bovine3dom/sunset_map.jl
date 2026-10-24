# Winter sunset map

Create an Arrow file of local sunset times for the world at H3 resolution 5.

Use the file with [H3-MON](https://github.com/bovine3dom/H3-MON) or whatever else u want

## Run

Use Julia 1.12 or later. Install the dependencies:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Create `sunset_2026_r5.arrow`:

```sh
julia --project=. --threads=auto sunset_map.jl
```

## Arrow columns (H3-MON schema)

| Column | Type | Meaning |
| --- | --- | --- |
| `index_lower` | `UInt32` | Lower 32 bits of the H3 index |
| `index_upper` | `UInt32` | Upper 32 bits of the H3 index |
| `value` | `Float32` | Local sunset time in minutes since midnight |

For example, `990.5` means 16:30:30 local time.
