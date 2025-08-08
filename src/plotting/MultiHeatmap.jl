"""
    plot_heatmap_axis(ax, data; kwargs...)

Plot a 2D heatmap on the given axis `ax` using PyPlot. Behaves like `plot_heatmap` but does not create a figure.

Keyword Arguments (all are optional):
- `x::Union{AbstractVector,Nothing}`: x-coordinates (default = 1:size(data,1))
- `y::Union{AbstractVector,Nothing}`: y-coordinates (default = 1:size(data,2))
- `cmap::String`: Colormap (default = "plasma")
- `xlabel::String`: x-axis label (default = "x")
- `ylabel::String`: y-axis label (default = "y")
- `label_fontsize::Real`: Font size for axis labels (default = 16)
- `tick_fontsize::Real`: Font size for tick labels (default = 14)
- `text_color::String`: Color for labels and title (default = "black")
- `tick_color::String`: Color for tick labels (default = "black")
- `show_xticks::Bool`: Show x-axis ticks (default = true)
- `show_yticks::Bool`: Show y-axis ticks (default = true)
- `show_xlabel::Bool`: Show x-axis label (default = true)
- `show_ylabel::Bool`: Show y-axis label (default = true)
- `xlim::Union{Tuple{<:Real,<:Real},Nothing}`: x-axis limits (default = nothing)
- `ylim::Union{Tuple{<:Real,<:Real},Nothing}`: y-axis limits (default = nothing)
- `title::Union{String,Nothing}`: Axis title (default = nothing)
- `title_fontsize::Real`: Font size for title (default = 18)
- `show_grid::Bool`: Show grid (default = false)
- `grid_alpha::Real`: Grid line transparency (default = 0.3)
- `grid_linewidth::Real`: Grid line width (default = 1.0)
- `show_colorbar::Bool`: Show colorbar (default = true)
- `colorbar_title::Union{String,Nothing}`: Colorbar label (default = nothing)
- `cb_label_fontsize::Real`: Colorbar label font size (default = 16)
- `cb_tick_fontsize::Real`: Colorbar tick font size (default = 14)
- `show_colorbar_ticks::Bool`: Show colorbar ticks (default = true)
- `show_colorbar_label::Bool`: Show colorbar label (default = true)
- `colorbar_orientation::String`: "vertical" or "horizontal" (default = "vertical")
- `colorbar_pad::Real`: Padding between plot and colorbar (default = 0.02)
- `colorbar_fraction::Real`: Fraction of original axes for colorbar (default = 0.046)
- `colorbar_shrink::Real`: Shrink factor for colorbar (default = 0.98)
- `colorbar_aspect::Union{Real,Nothing}`: Aspect ratio of colorbar (default = nothing)
- `colorbar_lim::Union{Tuple{<:Real,<:Real},Nothing}`: Set colorbar limits (vmin, vmax) (default = nothing)
- `colorbar_nticks::Union{Int,Nothing}`: Number of colorbar ticks (default = nothing)
- `overlay_lines::Vector{Tuple{Tuple{<:Real,<:Real}, Tuple{<:Real,<:Real}, String, String, Real}}`: Lines to overlay (default = [])
- `background_color::Union{String,Nothing}`: Background color for axis (default = nothing)

Returns:
- `im`: The AxesImage object
- `cb`: The Colorbar object (or `nothing` if not shown)
"""
function plot_heatmap_axis(ax::PyCall.PyObject, data::AbstractArray{<:Real,2}; kwargs...)
    # Size
    nx, ny = size(data)
    # Merge user options with defaults
    options = merge(Dict(
        :x => nothing,
        :y => nothing,
        :cmap => "plasma",
        :xlabel => "x",
        :ylabel => "y",
        :label_fontsize => 16,
        :tick_fontsize => 14,
        :text_color => "black",
        :tick_color => "black",
        :show_xticks => true,
        :show_yticks => true,
        :show_xlabel => true,
        :show_ylabel => true,
        :xlim => nothing,
        :ylim => nothing,
        :title => nothing,
        :title_fontsize => 18,
        :show_grid => false,
        :grid_alpha => 0.3,
        :grid_linewidth => 1.0,
        :show_colorbar => true,
        :colorbar_title => nothing,
        :cb_label_fontsize => 16,
        :cb_tick_fontsize => 14,
        :show_colorbar_ticks => true,
        :show_colorbar_label => true,
        :colorbar_orientation => "vertical",
        :colorbar_pad => 0.02,
        :colorbar_fraction => 0.046,
        :colorbar_shrink => 0.98,
        :colorbar_aspect => nothing,
        :colorbar_lim => nothing,
        :colorbar_nticks => nothing,
        :overlay_lines => [],
        :background_color => nothing
    ), Dict(kwargs))

    # Coordinate vectors
    x_vals = isnothing(options[:x]) ? (1:nx) : options[:x]
    y_vals = isnothing(options[:y]) ? (1:ny) : options[:y]

    # Figure reference for colorbar
    fig = ax.figure

    # Background color
    if !isnothing(options[:background_color])
        fig.patch.set_facecolor(options[:background_color])
        ax.set_facecolor(options[:background_color])
    end

    # Plot
    im = ax.imshow(data',
                   extent=(x_vals[1], x_vals[end], y_vals[1], y_vals[end]),
                   origin="lower", cmap=options[:cmap])

    # Labels
    if options[:show_xlabel]
        ax.set_xlabel(options[:xlabel], fontsize=options[:label_fontsize], color=options[:text_color])
    end
    if options[:show_ylabel]
        ax.set_ylabel(options[:ylabel], fontsize=options[:label_fontsize], color=options[:text_color])
    end

    # Ticks
    if !options[:show_xticks]
        ax.set_xticks([])
    else
        ax.tick_params(axis="x", labelsize=options[:tick_fontsize], colors=options[:tick_color])
    end
    if !options[:show_yticks]
        ax.set_yticks([])
    else
        ax.tick_params(axis="y", labelsize=options[:tick_fontsize], colors=options[:tick_color])
    end

    # Limits
    if !isnothing(options[:xlim])
        ax.set_xlim(options[:xlim])
    end
    if !isnothing(options[:ylim])
        ax.set_ylim(options[:ylim])
    end

    # Title
    if !isnothing(options[:title])
        ax.set_title(options[:title], fontsize=options[:title_fontsize], color=options[:text_color])
    end

    # Grid
    if options[:show_grid]
        ax.grid(true, alpha=options[:grid_alpha], linewidth=options[:grid_linewidth])
    end

    # Overlay lines
    for (p1, p2, color, style, width) in options[:overlay_lines]
        ax.plot([p1[1], p2[1]], [p1[2], p2[2]], color=color, linestyle=style, linewidth=width)
    end

    # Colorbar
    cb = nothing
    if options[:show_colorbar]
        cb_kwargs = Dict(
            :orientation => options[:colorbar_orientation],
            :pad => options[:colorbar_pad],
            :fraction => options[:colorbar_fraction],
            :shrink => options[:colorbar_shrink]
        )
        if !isnothing(options[:colorbar_aspect])
            cb_kwargs[:aspect] = options[:colorbar_aspect]
        end
        cb = fig.colorbar(im, ax=ax; cb_kwargs...)

        # Color limits
        if !isnothing(options[:colorbar_lim])
            vmin, vmax = options[:colorbar_lim]
            im.set_clim(vmin, vmax)
        end

        # Colorbar ticks
        if !options[:show_colorbar_ticks]
            cb.set_ticks([])
        elseif !isnothing(options[:colorbar_nticks])
            vmin, vmax = im.get_clim()
            ticks = range(vmin, vmax, length=options[:colorbar_nticks])
            cb.set_ticks(collect(ticks))
        end

        # Colorbar label
        if options[:show_colorbar_label]
            cb.set_label(options[:colorbar_title], fontsize=options[:cb_label_fontsize], color=options[:text_color])
        end
        cb.ax.tick_params(labelsize=options[:cb_tick_fontsize], colors=options[:tick_color])
    end

    return im, cb
end


"""
    process_axis_kwargs(kwargs, idx)

Extract and validate per-axis keyword arguments from a combined kwargs Dict.

# Arguments
- `kwargs::Dict`: Combined keyword arguments (any Dict type).
- `idx::Integer`: Index of the subplot (1-based) to extract.

# Returns
- `local_kwargs::Dict{Symbol,Any}`: Keyword arguments for this axis.

# Errors
- Throws an error if conflicting `key` and `key_all` are both present.
- Throws an error if any `key_all` value is not an AbstractVector or if its length < idx.
"""
function process_axis_kwargs(kwargs::Dict, idx::Integer)
    local_kwargs = Dict{Symbol,Any}()
    for (k,v) in kwargs
        s = String(k)
        if endswith(s, "_all")
            base = Symbol(chop(s, tail=4))  # remove `_all`
            if haskey(kwargs, base)
                error("Cannot specify both `$(base)` and `$(k)` in kwargs.")
            end
            if !(isa(v, AbstractVector))
                error("`$(k)` must be an AbstractVector, got $(typeof(v)).")
            end
            if length(v) < idx
                error("`$(k)` length ($(length(v))) is less than index $(idx).")
            end
            local_kwargs[base] = v[idx]
        elseif occursin("_all", s)
            # keys with `_all` in middle are not supported
            error("Keyword `$(k)` not recognized. Only suffix `_all` is supported for arrays.")
        else
            # non-array or non-_all key; keep for all
            local_kwargs[k] = v
        end
    end
    return local_kwargs
end

"""
    plot_heatmaps(data3d, nrows, ncols, indices; kwargs...)

Create a grid of heatmaps from a 3D dataset, with options for display, saving, and returning the figure.

# Arguments
- `data3d::AbstractArray{<:Real,3}`: 3D array of data to plot.
- `nrows::Integer`: Number of subplot rows.
- `ncols::Integer`: Number of subplot columns.
- `indices::AbstractVector{<:Integer}`: Indices of the third dimension to plot in each subplot, in row-major order.

# Keyword Arguments
All optional, passed via `kwargs...`. Supported kwargs:
- `:figsize`::Tuple — Figure size in inches, default `(6*ncols, 5*nrows)`.
- `:save`::Bool — Whether to save the figure, default `false`.
- `:figname`::String — Filename to save if `save=true`, default `"heatmaps.png"`.
- `:dpi`::Int — Resolution for saving, default `300`.
- `:bbox_inches`::String — Bounding box for saving, default `"tight"`.
- `:transparent`::Bool — Save with transparent background, default `false`.
- `:return_objects`::Bool — Return `(fig, axes)` instead of displaying, default `false`.
- For any `param`, you may specify `:param => value` to apply to all subplots, or `:param_all => vector` (length ≥ number of subplots) for per-axis values.

Supported `param` keys mirror those of `plot_heatmap_axis`, e.g. `:title`, `:xlabel`, `:cmap`, etc.

# Returns
- If `return_objects=true`, returns `(fig, axes)`.
- Otherwise, displays or saves the figure as specified.
"""
function plot_heatmaps(
    data3d::AbstractArray{<:Real,3},
    nrows::Integer,
    ncols::Integer,
    indices::AbstractVector{<:Integer};
    kwargs...
)
    # Validate number of subplots
    nplots = nrows * ncols
    nin = length(indices)
    if nin > nplots
        error("Number of indices ($(nin)) exceeds number of subplots ($(nplots)).")
    elseif nin < nplots
        @warn "Number of indices ($(nin)) is less than number of subplots ($(nplots)). Remaining subplots will be empty."
    end

    # Merge default options with user kwargs
    defaults = Dict(
        :figsize => (6*ncols, 5*nrows),
        :save => false,
        :figname => "heatmaps.png",
        :dpi => 300,
        :bbox_inches => "tight",
        :transparent => false,
        :return_objects => false
    )
    kw = merge(defaults, Dict(kwargs...))

    # Create figure and axes using merged figsize
    fig, axarr = subplots(nrows, ncols; figsize=kw[:figsize])
  if nrows ==1 || ncols == 1
        axes = vec(axarr)
    else
        axes = [axarr[r, c] for r in 1:nrows for c in 1:ncols]
    end

    # Plot each requested slice
    for i in 1:min(nin, nplots)
        idx = indices[i]
        ax = axes[i]
        data2d = data3d[:, :, idx]
        local_kwargs = process_axis_kwargs(kw, i)
        plot_heatmap_axis(ax, data2d; local_kwargs...)
    end

    # Hide unused axes
    if nin < nplots
        for i in (nin+1):nplots
            axes[i].axis("off")
        end
    end

    # Save if requested using merged options
    if kw[:save]
        fig.savefig(kw[:figname]; dpi=kw[:dpi], bbox_inches=kw[:bbox_inches], transparent=kw[:transparent])
        println("Figure saved as $(kw[:figname])")
    end

    # Return or display
    if kw[:return_objects]
        return fig, axes
    elseif !kw[:save]
        display(fig)
    end
end


"""
    plot_heatmaps(data4d, nrows, ncols, indices, timeslot; kwargs...)

Create a grid of heatmaps from a 4D dataset with variable index and timeslot.

# Arguments
- `data4d::AbstractArray{<:Real,4}`: 4D data (x, y, variable, time).
- `nrows::Integer`, `ncols::Integer`: Subplot grid dimensions.
- `indices::Union{Integer, AbstractVector{<:Integer}}`: Variable indices.
- `timeslot::Union{Integer, AbstractVector{<:Integer}}`: Time indices.

# Keyword Arguments
Same as 3D version (`:figsize`, `:save`, `:figname`, `:dpi`, `:bbox_inches`, `:transparent`, `:return_objects`) and any `param` or `:param_all` keys for per-axis options.

# Behavior
- Normalizes `indices` and `timeslot` to vectors of equal length (expanding scalars).
- Validates number of subplot requests against `nrows*ncols`.
- Uses each `(idx, tslot)` pair to slice `data4d[:, :, idx, tslot]`.
- Delegates plotting to `plot_heatmap_axis` with processed kwargs.

# Returns
- Same as 3D version.
"""
function plot_heatmaps(
    data4d::AbstractArray{<:Real,4},
    nrows::Integer,
    ncols::Integer,
    indices::Union{Integer, AbstractVector{<:Integer}},
    timeslot::Union{Integer, AbstractVector{<:Integer}};
    kwargs...
)
    # Total slots and normalize inputs
    nplots = nrows * ncols
    idx_vec = isa(indices, Integer) ? [indices] : collect(indices)
    tslot_vec = isa(timeslot, Integer) ? [timeslot] : collect(timeslot)

    if length(idx_vec) == 1 && length(tslot_vec) > 1
        idx_vec = fill(idx_vec[1], length(tslot_vec))
    elseif length(tslot_vec) == 1 && length(idx_vec) > 1
        tslot_vec = fill(tslot_vec[1], length(idx_vec))
    end
    nreq = length(idx_vec)

    if nreq > nplots
        error("Number of subplot requests ($(nreq)) exceeds available subplots ($(nplots)).")
    elseif nreq < nplots
        @warn "Number of subplot requests ($(nreq)) is less than available subplots ($(nplots)). Remaining subplots will be empty."
    end

    # Merge defaults
    defaults = Dict(
        :figsize => (6*ncols, 5*nrows),
        :save => false,
        :figname => "heatmaps.png",
        :dpi => 300,
        :bbox_inches => "tight",
        :transparent => false,
        :return_objects => false
    )
    kw = merge(defaults, Dict(kwargs...))

    # Create figure
    fig, axarr = subplots(nrows, ncols; figsize=kw[:figsize])
    if nrows ==1 || ncols == 1
        axes = vec(axarr)
    else
        axes = [axarr[r, c] for r in 1:nrows for c in 1:ncols]
    end
    # Plot requested slices
    for i in 1:min(nreq, nplots)
        ax = axes[i]
        data2d = data4d[:, :, idx_vec[i], tslot_vec[i]]
        local_kwargs = process_axis_kwargs(kw, i)
        plot_heatmap_axis(ax, data2d; local_kwargs...)
    end

    # Hide extras
    if nreq < nplots
        for i in (nreq+1):nplots
            axes[i].axis("off")
        end
    end

    # Save and return/display
    if kw[:save]
        fig.savefig(kw[:figname]; dpi=kw[:dpi], bbox_inches=kw[:bbox_inches], transparent=kw[:transparent])
        println("Figure saved as $(kw[:figname])")
    end
    if kw[:return_objects]
        return fig, axes
    elseif !kw[:save]
        display(fig)
    end
end