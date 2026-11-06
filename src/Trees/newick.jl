"""
    newick(nodes; label = nothing, annotate = nothing)
    newick(io::IO, nodes; label = nothing, annotate = nothing)

Returns the nodes of the node collection `nodes` in Newick format. Nodes whose parent is not in `nodes` are treated as roots. Branch lengths are lifetimes. Requires a sorted tree.

Nodes are labelled by their ids, or by `label(node)` if given. If given, `annotate(node)` returns a named tuple (or `nothing`) to be added as an annotation of the form `[&key=value,...]`.
"""
function newick(io::IO, nodes; label = nothing, annotate = nothing)
    ensure_sorted(nodes)
    tree = get_tree(nodes)
    keep = falses(length(tree))
    keep[collect(keys(nodes))] .= true

    for id in findall(keep)
        p = parent(tree, id)
        p == 0 || !keep[p] || continue

        # Iterative, as lineages can be very deep. `i > 0` opens node `i`, `-i` closes it, `0` separates siblings.
        stack = [ id ]
        while !isempty(stack)
            i = pop!(stack)
            cs = i > 0 ? filter(c -> keep[c], child_ids(tree, i)) : Int[]
            if i == 0
                print(io, ',')
            elseif !isempty(cs)
                print(io, '(')
                push!(stack, -i)
                foreach(c -> push!(stack, c, 0), reverse(cs))
                pop!(stack)
            else
                i < 0 && print(io, ')')
                write_node(io, abs(i), tree[abs(i)], label, annotate)
            end
        end

        println(io, ';')
    end
end

newick(nodes; kwargs...) = sprint(io -> newick(io, nodes; kwargs...))

function write_node(io, id, node, label, annotate)
    # Labels with spaces or special characters are quoted
    s = string(something(isnothing(label) ? id : label(node), ""))
    print(io, occursin(r"[\s()\[\]':;,]", s) ? "'" * replace(s, "'" => "''") * "'" : s)

    ann = isnothing(annotate) ? nothing : annotate(node)
    if !isnothing(ann) && !isempty(ann)
        print(io, "[&")
        join(io, (string(k, '=', v isa Real ? v : "\"$v\"") for (k, v) in pairs(ann)), ',')
        print(io, ']')
    end

    print(io, ':', lifetime(node))
end
