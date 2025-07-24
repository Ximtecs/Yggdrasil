"""
    plot_line_axis(ax, data; kwargs...)

Plot a 1D line on the given axis `ax` using PyPlot.

Keyword Arguments (all optional):
- `x::Union{AbstractVector{<:Real},Nothing}`: x-values (default = `1:length(data)`).
- `color::String`: Line color (default = "black").
- `linewidth::Real`: Line thickness (default = 2.0).
- `linestyle::String`: Line style (default = "-").
- `marker::Union{String,Nothing}`: Marker style (default = `nothing`).
- `markersize::Real`: Marker size (default = 6).
- `markevery::Union{Int, AbstractVector{<:Integer}, Nothing}`: Marker frequency (default = `nothing`).
- `label::Union{String,Nothing}`: Legend label for the line (default = `nothing`).

Plots the data vector against provided x-values on `ax`.
"""
function plot_line_axis(ax::PyCall.PyObject, data::AbstractVector{<:Real}; kwargs...)
    # Default line and x-options
    defaults = Dict(
        :x => nothing,
        :color => "black",
        :linewidth => 2.0,
        :linestyle => "-",
        :marker => nothing,
        :markersize => 6,
        :markevery => nothing,
        :label => nothing
    )
    opts = merge(defaults, Dict(kwargs...))

    # Determine x-coordinates
    x_vals = isnothing(opts[:x]) ? collect(1:length(data)) : opts[:x]

    # Plot the line
    ax.plot(
        x_vals,
        data;
        color=opts[:color],
        linewidth=opts[:linewidth],
        linestyle=opts[:linestyle],
        marker=opts[:marker],
        markersize=opts[:markersize],
        markevery=opts[:markevery],
        label=opts[:label]
    )
end

"""
    plot_line(data; kwargs...)

Create a simple line plot for a 1D data vector.

Arguments:
- `data::AbstractVector{<:Real}`: Data to plot (y-values).

Keyword Arguments:
- `figsize::Tuple{<:Real,<:Real}`: Figure size in inches (default = (8,4)).
- `facecolor::String`: Figure background color (default = "white").
- `ax_facecolor::String`: Axis background color (default = "white").
- `save::Bool`: Whether to save the figure (default = false).
- `figname::String`: Filename to save if `save=true` (default = "lineplot.png").
- `dpi::Int`: DPI for saving (default = 300).
- `bbox_inches::String`: Bounding box for saving (default = "tight").
- `transparent::Bool`: Save with transparent background (default = false).
- `return_objects::Bool`: Return `(fig, ax)` instead of display (default = false).
- `xlabel::String`: X-axis label (default = "x").
- `ylabel::String`: Y-axis label (default = "y").
- `title::Union{String,Nothing}`: Plot title (default = nothing).
- `grid::Bool`: Show grid (default = true).
- `grid_color::String`: Grid line color (default = "grey").
- `grid_linestyle::String`: Grid line style (default = "--").
- `grid_linewidth::Real`: Grid line width (default = 1.0).
- `xlabel_fontsize::Real`: X-axis label font size (default = 12).
- `ylabel_fontsize::Real`: Y-axis label font size (default = 12).
- `title_fontsize::Real`: Title font size (default = 14).
- `tick_fontsize::Real`: Tick label font size (default = 10).
- `legend::Bool`: Show legend if any labels provided (default = false).
- `legend_loc::String`: Legend position (default = "best").
- `legend_fontsize::Real`: Legend font size (default = 10).
- Any other kwargs are forwarded to `plot_line_axis`, including `:x` and `:label`.
"""
function plot_line(data::AbstractVector{<:Real}; kwargs...)
    # Merge default figure and axis options
    defaults = Dict(
        :figsize => (8,4),
        :facecolor => "white",
        :ax_facecolor => "white",
        :save => false,
        :figname => "lineplot.png",
        :dpi => 300,
        :bbox_inches => "tight",
        :transparent => false,
        :return_objects => false,
        :xlabel => "x",
        :ylabel => "y",
        :title => nothing,
        :grid => true,
        :grid_color => "grey",
        :grid_linestyle => "--",
        :grid_linewidth => 1.0,
        :xlabel_fontsize => 12,
        :ylabel_fontsize => 12,
        :title_fontsize => 14,
        :tick_fontsize => 10,
        :legend => false,
        :legend_loc => "best",
        :legend_fontsize => 10
    )
    kw = merge(defaults, Dict(kwargs...))

    # Create figure and axis
    fig, ax = subplots(figsize=kw[:figsize])
    fig.patch.set_facecolor(kw[:facecolor])
    ax.set_facecolor(kw[:ax_facecolor])

    # Plot on axis (forwards all kwargs, including x and label)
    plot_line_axis(ax, data; kwargs...)

    # Apply axis labels, title, and font sizes
    ax.set_xlabel(kw[:xlabel], fontsize=kw[:xlabel_fontsize])
    ax.set_ylabel(kw[:ylabel], fontsize=kw[:ylabel_fontsize])
    if kw[:title] !== nothing
        ax.set_title(kw[:title], fontsize=kw[:title_fontsize])
    end
    if kw[:grid]
        ax.grid(true, color=kw[:grid_color], linestyle=kw[:grid_linestyle], linewidth=kw[:grid_linewidth])
    end

    # Tick labels font size
    ax.tick_params(labelsize=kw[:tick_fontsize])

    # Legend
    if kw[:legend]
        ax.legend(loc=kw[:legend_loc], fontsize=kw[:legend_fontsize])
    end

    # Save if requested
    if kw[:save]
        fig.savefig(kw[:figname]; dpi=kw[:dpi], bbox_inches=kw[:bbox_inches], transparent=kw[:transparent])
        println("Figure saved as $(kw[:figname])")
    end

    # Return or display
    if kw[:return_objects]
        return fig, ax
    elseif !kw[:save]
        display(fig)
    end
end


#------------------------------------------------------------------------------------------------

"""
    plot_line(data; kwargs...)

Plot multiple lines from either a 2D data matrix or a vector of line vectors.

# Arguments
- `data::Union{AbstractMatrix{<:Real}, AbstractVector{<:AbstractVector{<:Real}}}`: 2D matrix with columns as lines, or a vector of 1D vectors.

# Keyword Arguments
- `indices::Union{Nothing, AbstractVector{<:Integer}}`: Lines to plot (default = `nothing`, plots all).
- Figure/axis-level options (see `plot_line` defaults):
  `:figsize, :facecolor, :ax_facecolor, :save, :figname, :dpi, :bbox_inches, :transparent, :return_objects,`
  `:xlabel, :ylabel, :title, :grid, :grid_color, :grid_linestyle, :grid_linewidth,`
  `:xlabel_fontsize, :ylabel_fontsize, :title_fontsize, :tick_fontsize, :legend, :legend_loc, :legend_fontsize`.
- Line-specific styling via `_all` suffix or uniform: `:x, :color, :linewidth, :linestyle, :marker, :markersize, :markevery, :label`.

# Behavior
1. Determine `nlines` and accessor `get_y(i)`. Errors otherwise.
2. Normalize `indices`: `nothing`→all lines, else must be vector of Int.
3. Validate all indices in `1:nlines`.
4. Inject defaults:
   - If no `:color_all`, set it to first `length(idxs)` of `default_color_palette`, cycling.
   - If no `:label_all`, set `"line 1"`, …, `"line N"`.
5. Create figure/axis.
6. For each selected line `i`:
   - `local_kwargs = process_axis_kwargs(kw, i)`
   - `plot_line_axis(ax, get_y(idx); local_kwargs...)
7. Apply figure-level styling, grid, ticks, legend.
8. Save or return/display.

# Returns
- `(fig, ax)` if `:return_objects=true`, else shows/saves.
"""
function plot_line(data::Union{AbstractMatrix{<:Real}, AbstractVector{<:AbstractVector{<:Real}}}; kwargs...)
    # Default options
    defaults = Dict(
        :indices => nothing,
        :figsize => (8,4),
        :facecolor => "white",
        :ax_facecolor => "white",
        :save => false,
        :figname => "lineplot.png",
        :dpi => 300,
        :bbox_inches => "tight",
        :transparent => false,
        :return_objects => false,
        :xlabel => "x",
        :ylabel => "y",
        :title => nothing,
        :grid => true,
        :grid_color => "grey",
        :grid_linestyle => "--",
        :grid_linewidth => 1.0,
        :xlabel_fontsize => 12,
        :ylabel_fontsize => 12,
        :title_fontsize => 14,
        :tick_fontsize => 10,
        :legend => true,
        :legend_loc => "best",
        :legend_fontsize => 10
    )
    kw = merge(defaults, Dict(kwargs...))

    # Determine nlines and y-accessor
    nlines = 0
    get_y = i -> error()
    if isa(data, AbstractMatrix)
        nlines = size(data, 2)
        get_y = i -> data[:, i]
    elseif isa(data, AbstractVector)
        nlines = length(data)
        get_y = i -> data[i]
    else
        error("`data` must be a matrix or a vector of vectors.")
    end

    # Indices
    idxs = kw[:indices] === nothing ? collect(1:nlines) : kw[:indices]
    if !(kw[:indices] === nothing || isa(idxs, AbstractVector{<:Integer}))
        error("`indices` must be a Vector{Int} or nothing.")
    end
    for c in idxs
        if c < 1 || c > nlines
            error("Index $(c) out of bounds (1:$(nlines)).")
        end
    end

    # Inject default color_all
    if !haskey(kw, :color_all)
        kw[:color_all] = [ default_color_palette[mod1(i, length(default_color_palette))] for i in 1:length(idxs) ]
    end
    # Inject default label_all
    if !haskey(kw, :label_all)
        kw[:label_all] = ["line " * string(i) for i in 1:length(idxs)]
    end

    # Create figure/axis
    fig, ax = subplots(figsize=kw[:figsize])
    fig.patch.set_facecolor(kw[:facecolor])
    ax.set_facecolor(kw[:ax_facecolor])

    # Plot lines
    for (i, idx) in enumerate(idxs)
        local_kwargs = process_axis_kwargs(kw, i)
        plot_line_axis(ax, get_y(idx); local_kwargs...)
    end

    # Figure-level styling
    ax.set_xlabel(kw[:xlabel], fontsize=kw[:xlabel_fontsize])
    ax.set_ylabel(kw[:ylabel], fontsize=kw[:ylabel_fontsize])
    if kw[:title] !== nothing
        ax.set_title(kw[:title], fontsize=kw[:title_fontsize])
    end
    if kw[:grid]
        ax.grid(true, color=kw[:grid_color], linestyle=kw[:grid_linestyle], linewidth=kw[:grid_linewidth])
    end
    ax.tick_params(labelsize=kw[:tick_fontsize])
    if kw[:legend]
        ax.legend(loc=kw[:legend_loc], fontsize=kw[:legend_fontsize])
    end

    # Save or return/display
    if kw[:save]
        fig.savefig(kw[:figname]; dpi=kw[:dpi], bbox_inches=kw[:bbox_inches], transparent=kw[:transparent])
    end
    if kw[:return_objects]
        return fig, ax
    elseif !kw[:save]
        display(fig)
    end
end
