"""
    plot_heatmap(data; kwargs...)

Plot a 2D heatmap using PyPlot. Only `data` is required; all other options are set via keyword arguments.

# Keyword Arguments

## Data layout
- `x::Vector`: x-coordinates (default = 1:size(data,1))
- `y::Vector`: y-coordinates (default = 1:size(data,2))
- `cmap::String`: Colormap (default = "plasma")

## Axis labelling and ticks
- `xlabel::String`: x-axis label (default = "x")
- `ylabel::String`: y-axis label (default = "y")
- `label_fontsize::Int`: Axis label font size (default = 16)
- `tick_fontsize::Int`: Tick label font size (default = 14)
- `text_color::String`: Colour of all labels and title (default = "black")
- `tick_color::String`: Colour of axis and colourbar tick labels (default = "black")
- `show_xticks::Bool`: Show x-axis ticks (default = true)
- `show_yticks::Bool`: Show y-axis ticks (default = true)
- `show_xlabel::Bool`: Show x-axis label (default = true)
- `show_ylabel::Bool`: Show y-axis label (default = true)
- `xlim::Union{Nothing, Tuple}`: Set x-axis limits (default = nothing)
- `ylim::Union{Nothing, Tuple}`: Set y-axis limits (default = nothing)

## Title
- `title::Union{String,Nothing}`: Optional plot title (default = nothing)
- `title_fontsize::Int`: Title font size (default = 18)

## Grid
- `show_grid::Bool`: Show grid aligned with ticks (default = false)
- `grid_alpha::Float64`: Transparency of grid lines (default = 0.3)
- `grid_linewidth::Float64`: Width of grid lines (default = 1.0)

## Colourbar
- `show_colorbar::Bool`: Show colourbar (default = true)
- `colorbar_title::String`: Label for the colourbar (default = "Field")
- `cb_label_fontsize::Int`: Colourbar label font size (default = 16)
- `cb_tick_fontsize::Int`: Colourbar tick font size (default = 14)
- `show_colorbar_ticks::Bool`: Show colourbar ticks (default = true)
- `show_colorbar_label::Bool`: Show colourbar label (default = true)
- `colorbar_orientation::String`: "vertical" or "horizontal" (default = "vertical")
- `colorbar_pad::Float64`: Padding between plot and colourbar (default = 0.02)
- `colorbar_fraction::Float64`: Fraction of original axes for colourbar (default = 0.046)
- `colorbar_shrink::Float64`: Shrink factor for colourbar (default = 0.9)
- `colorbar_aspect::Union{Nothing, Real}`: Aspect ratio of colourbar (default = nothing)
- `colorbar_lim::Union{Nothing, Tuple{<:Real,<:Real}}`: Set colourbar limits as (vmin, vmax) (default = nothing)
- `colorbar_nticks::Union{Nothing, Int}`: Number of colourbar ticks (default = nothing)

## Overlay
- `overlay_lines`: Vector of line specs (p1, p2, color, linestyle, linewidth)

## Layout and saving
- `figsize::Tuple`: Figure size (default = (8, 6))
- `dpi::Int`: Resolution for saving (default = 300)
- `bbox_inches::String`: Bounding box for saving (default = "tight")
- `save::Bool`: Whether to save the figure (default = false)
- `figname::String`: Filename if `save=true` (default = "field_plot.png")
- `transparent::Bool`: Save with transparent background (default = false)
- `background_color::Union{String, Nothing}`: Background colour for figure and axes (default = nothing)
- `return_objects::Bool`: Return `fig`, `ax`, `im` instead of displaying (default = false)
"""
function plot_heatmap(data::AbstractArray{<:Real,2}; kwargs...)

    nx, ny = size(data)
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
        :figsize => (8, 6),
        :dpi => 300,
        :bbox_inches => "tight",
        :save => false,
        :figname => "field_plot.png",
        :transparent => false,
        :background_color => nothing,
        :return_objects => false
    ), Dict(kwargs))

    x_vals = isnothing(options[:x]) ? (1:nx) : options[:x]
    y_vals = isnothing(options[:y]) ? (1:ny) : options[:y]

    fig, ax = subplots(figsize=options[:figsize])

    if !isnothing(options[:background_color])
        fig.patch.set_facecolor(options[:background_color])
        ax.set_facecolor(options[:background_color])
    end

    im = ax.imshow(data',
                   extent=(x_vals[1], x_vals[end], y_vals[1], y_vals[end]),
                   origin="lower", cmap=options[:cmap])

    if options[:show_xlabel]
        ax.set_xlabel(options[:xlabel], fontsize=options[:label_fontsize], color=options[:text_color])
    end
    if options[:show_ylabel]
        ax.set_ylabel(options[:ylabel], fontsize=options[:label_fontsize], color=options[:text_color])
    end
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

    if !isnothing(options[:xlim])
        ax.set_xlim(options[:xlim])
    end
    if !isnothing(options[:ylim])
        ax.set_ylim(options[:ylim])
    end

    if !isnothing(options[:title])
        ax.set_title(options[:title], fontsize=options[:title_fontsize], color=options[:text_color])
    end

    if options[:show_grid]
        ax.grid(true, alpha=options[:grid_alpha], linewidth=options[:grid_linewidth])
    end

    # Overlay lines
    for (p1, p2, color, style, width) in options[:overlay_lines]
        ax.plot([p1[1], p2[1]], [p1[2], p2[2]], color=color, linestyle=style, linewidth=width)
    end

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

        if !isnothing(options[:colorbar_lim])
            vmin, vmax = options[:colorbar_lim]
            im.set_clim(vmin, vmax)
        else
            vmin, vmax = im.get_clim()
        end

        if !options[:show_colorbar_ticks]
            cb.set_ticks([])
        elseif !isnothing(options[:colorbar_nticks])
            ticks = range(vmin, vmax, length=options[:colorbar_nticks])
            cb.set_ticks(collect(ticks))
        end

        if options[:show_colorbar_label]
            cb.set_label(options[:colorbar_title], fontsize=options[:cb_label_fontsize], color=options[:text_color])
        end
        cb.ax.tick_params(labelsize=options[:cb_tick_fontsize], colors=options[:tick_color])
    end

    if options[:save]
        fig.savefig(options[:figname], dpi=options[:dpi], bbox_inches=options[:bbox_inches], transparent=options[:transparent])
        println("Figure saved as $(options[:figname])")
    end

    if options[:return_objects]
        return fig, ax, im, options[:show_colorbar] ? cb : nothing
    end

    if !options[:save] && !options[:return_objects]
        display(fig)
    end

end


function plot_heatmap(data::AbstractArray{<:Real, 3}, var::Int; kwargs...)
    reduced = drop_unit_dims(data)
    if ndims(reduced) != 3
        if ndims(reduced) == 2
            return plot_heatmap(reduced; kwargs...)
        else
            error("Invalid data dimensions: $(ndims(reduced)). Expected 2, 3")
        end
    end

    options = merge(Dict(:var_idx => 3), Dict(kwargs))
    vi = options[:var_idx]

    if vi < 1 || vi > 3
        error("Invalid var_idx: $(vi). Must be 1, 2, or 3.")
    end

    idx = Vector{Any}(undef, 3); fill!(idx, Colon()) 
    idx[vi] = var

    slice = view(reduced, idx...)
    return plot_heatmap(slice; kwargs...)
end

function plot_heatmap(data::AbstractArray{<:Real, 3}; kwargs...)
    reduced = drop_unit_dims(data)
    if ndims(reduced) != 3
        return plot_heatmap(reduced; kwargs...)
    else
        error("3D data must either have collapsable dimension (dim with single index) or have the variable spcified: plot_heatmap(data, var)")
    end
end



function plot_heatmap(data::AbstractArray{<:Real, 4}, var::Int, t::Int; kwargs...)
    reduced = drop_unit_dims(data)
    if ndims(reduced) != 4
        if ndims(reduced) == 3
            return plot_heatmap(reduced, var; kwargs...)
        elseif ndims(reduced) == 2
            return plot_heatmap(reduced; kwargs...)
        else
            error("Invalid data dimensions: $(ndims(reduced)). Expected 2, 3, or 4.")
        end
    end

    options = merge(Dict(:var_idx => 3, :time_idx => 4), Dict(kwargs))
    vi = options[:var_idx]
    ti = options[:time_idx]

    if vi < 1 || vi > 4 || ti < 1 || ti > 4 || vi == ti
        error("Invalid var_idx ($(vi)) or time_idx ($(ti)) for 4D array.")
    end

    idx = Vector{Any}(undef, 4); fill!(idx, Colon()) 
    idx[vi] = var
    idx[ti] = t

    slice = view(reduced, idx...)
    return plot_heatmap(slice; kwargs...)
end


function plot_heatmap(data::AbstractArray{<:Real, 4}, var::Int; kwargs...)
    reduced = drop_unit_dims(data)
    if ndims(reduced) != 4
        if ndims(reduced) == 3
            return plot_heatmap(reduced, var; kwargs...)
        elseif ndims(reduced) == 2
            return plot_heatmap(reduced; kwargs...)
        else
            error("Invalid data dimensions: $(ndims(reduced)). Expected 2, 3")
        end
    else
        error("4D data must either have collapsable time dimension or have both variable and time: plot_heatmap(data, var, t)")
    end
end

function plot_heatmap(data::AbstractArray{<:Real, 4}; kwargs...)
    reduced = drop_unit_dims(data)
    if ndims(reduced) != 4
        if ndims(reduced) == 2
            return plot_heatmap(reduced; kwargs...)
        else
            error("Invalid data dimensions: $(ndims(reduced)). Expected 2.")
        end
    else
        error("4D data must either have collapsable time dimension or have both variable and time: plot_heatmap(data, var, t)")
    end
end