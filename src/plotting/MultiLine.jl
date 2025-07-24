"""
    plot_lines(data, nrows, ncols; kwargs...)

Plot multiple lines into a grid of subplots, assigning each line to a specific subplot.

Arguments:
- `data::Union{AbstractMatrix{<:Real}, AbstractVector{<:AbstractVector{<:Real}}}`: 2D matrix (columns are lines) or a vector of 1D real vectors.
- `nrows::Int`, `ncols::Int`: Number of rows and columns in the subplot grid.

Keyword Arguments:
- `indices::Union{Nothing, Vector{Int}}`: Lines to plot (default `nothing` ⇒ all lines).
- `axis_map::Union{Nothing, Vector{Int}}`: Subplot index (1 to `nrows*ncols`) for each line; if `nothing`, requires one line per subplot and maps in row-major order.
- Figure-level options (as in `plot_line`):
  `:figsize`, `:facecolor`, `:ax_facecolor`, `:save`, `:figname`, `:dpi`, `:bbox_inches`, `:transparent`,
  `:return_objects`, `:xlabel`, `:ylabel`, `:title`, `:grid`, `:grid_color`, `:grid_linestyle`, `:grid_linewidth`,
  `:xlabel_fontsize`, `:ylabel_fontsize`, `:title_fontsize`, `:tick_fontsize`, `:legend`, `:legend_loc`, `:legend_fontsize`.
- Line-specific styling forwarded via `process_axis_kwargs`: `_all` suffix or uniform for `:x`, `:color`, `:linewidth`, `:linestyle`, `:marker`, `:markersize`, `:markevery`, `:label`.

Returns:
- `(fig, axes)` if `:return_objects=true`, otherwise displays (and saves if `:save=true`).
"""
function plot_lines(data::Union{AbstractMatrix{<:Real}, AbstractVector{<:AbstractVector{<:Real}}},
                    nrows::Int, ncols::Int; kwargs...)
    # Merge defaults
    defaults = Dict(
        :indices           => nothing,
        :axis_map          => nothing,
        :figsize           => (8*ncols, 5*nrows),
        :facecolor         => "white",
        :ax_facecolor      => "white",
        :save              => false,
        :figname           => "lineplots.png",
        :dpi               => 300,
        :bbox_inches       => "tight",
        :transparent       => false,
        :return_objects    => false,
        :xlabel            => "x",
        :ylabel            => "y",
        :fig_title         => nothing,
        :sub_titles        => nothing,
        :grid              => true,
        :grid_color        => "grey",
        :grid_linestyle    => "--",
        :grid_linewidth    => 1.0,
        :xlabel_fontsize   => 12,
        :ylabel_fontsize   => 12,
        :subplot_title_fontsize => 14,
        :suptitle_fontsize => 16,
        :tick_fontsize     => 10,
        :legend            => true,
        :legend_loc        => "best",
        :legend_fontsize   => 10
    )
    kw = merge(defaults, Dict(kwargs...))

    # Determine line count and accessor
    nlines = isa(data,AbstractMatrix) ? size(data,2) : length(data)
    get_y = i -> isa(data,AbstractMatrix) ? data[:,i] : data[i]

    # Normalize and validate indices
    idxs = kw[:indices]===nothing ? collect(1:nlines) : kw[:indices]
    if !(kw[:indices]===nothing || isa(idxs,Vector{Int}))
        error("`indices` must be a Vector{Int} or nothing.")
    end
    for c in idxs
        if c<1||c>nlines
            error("Index $(c) out of bounds (1:$(nlines)).")
        end
    end

    # Build or validate axis_map
    nplots = nrows*ncols
    amap = kw[:axis_map]
    if amap===nothing
        if length(idxs)!=nplots
            error("With no axis_map, number of lines ($(length(idxs))) must equal subplots ($(nplots)).")
        end
        amap = collect(1:nplots)
    else
        if !(isa(amap,Vector{Int}) && length(amap)==length(idxs))
            error("`axis_map` must be Vector{Int} matching number of lines ($(length(idxs))).")
        end
        for p in amap
            if p<1||p>nplots
                error("axis_map entry $(p) out of range (1:$(nplots)).")
            end
        end
    end

    # Inject default color_all per subplot
    if !haskey(kw,:color_all)
        color_all = Vector{String}(undef, length(idxs))
        for sub in 1:nplots
            lines_in_sub = findall(x->x==sub, amap)
            for (pos,j) in enumerate(lines_in_sub)
                color_all[j] = default_color_palette[mod1(pos,length(default_color_palette))]
            end
        end
        kw[:color_all] = color_all
    end
    # Inject default labels per subplot
    if !haskey(kw,:label_all)
        kw[:label_all] = ["line "*string(i) for i in 1:length(idxs)]
    end

    # Create subplots grid
    fig, axarr = subplots(nrows,ncols; figsize=kw[:figsize])
    axes = [axarr[r,c] for r in 1:nrows for c in 1:ncols]
    fig.patch.set_facecolor(kw[:facecolor])
    for ax in axes
        ax.set_facecolor(kw[:ax_facecolor])
    end


   # Subplot titles
    if kw[:sub_titles] !== nothing
        if !(isa(kw[:sub_titles], Vector) && length(kw[:sub_titles]) == nplots)
            println(isa(kw[:sub_titles], Vector), " ---- ", length(kw[:sub_titles]) == nplots)
            error("`sub_titles` must be a Vector with length equal to the number of subplots ($(nplots)).")
        end
        for p in 1:nplots
            axes[p].set_title(kw[:sub_titles][p], fontsize=kw[:subplot_title_fontsize])
        end
    end



    # Plot each line in its assigned subplot
    for (i, idx) in enumerate(idxs)
        ax = axes[amap[i]]
        local_kw = process_axis_kwargs(kw, i)
        plot_line_axis(ax, get_y(idx); local_kw...)
    end

    # Overall figure title
    if kw[:fig_title] !== nothing
        fig.suptitle(kw[:fig_title], fontsize=kw[:suptitle_fontsize])
    end

    # Apply styling to all subplots
    for ax in axes
        ax.set_xlabel(kw[:xlabel], fontsize=kw[:xlabel_fontsize])
        ax.set_ylabel(kw[:ylabel], fontsize=kw[:ylabel_fontsize])
        if kw[:grid]
            ax.grid(true, color=kw[:grid_color], linestyle=kw[:grid_linestyle], linewidth=kw[:grid_linewidth])
        end
        ax.tick_params(labelsize=kw[:tick_fontsize])
        if kw[:legend]
            ax.legend(loc=kw[:legend_loc], fontsize=kw[:legend_fontsize])
        end
    end

    # Save or return/display
    if kw[:save]
        fig.savefig(kw[:figname]; dpi=kw[:dpi], bbox_inches=kw[:bbox_inches], transparent=kw[:transparent])
    end
    if kw[:return_objects]
        return fig, axes
    elseif !kw[:save]
        display(fig)
    end
end
