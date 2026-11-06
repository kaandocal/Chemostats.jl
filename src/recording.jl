# Recording into `chem.tree` (a `Trees.Tree`, or `nothing` if not recording).
# Calls during a simulation must hold `int.tree_lock`.

record_init!(::Nothing, cells) = nothing
function record_init!(tree::Tree, cells)
    for cell in cells
        Trees.cell_init!(tree, (cell,), get_curr_t(cell))
    end
end

record_division!(::Nothing, parent, children, t, data) = nothing
record_division!(tree::Tree, parent, children, t, data) = Trees.cell_divide!(tree, parent, children, t, data)

record_end!(::Nothing, cell, t, data) = nothing
record_end!(tree::Tree, cell, t, data) = Trees.cell_die!(tree, cell, t, data)

# Ends `cell` at its current time if it is still living in the tree
record_drop!(::Nothing, cell) = nothing
function record_drop!(tree::Tree, cell)
    haskey(tree.living, cell) && Trees.cell_die!(tree, cell, get_curr_t(cell), get_data(cell))
    nothing
end

# Keeps the current population when pruning, so that snapshots at this time stay exact
mark_tree!(chem) = (isnothing(chem.tree) || Trees.mark!(chem.tree); nothing)

# An extinct population would be pruned away entirely, so it is kept
function prune_tree!(chem)
    (; tree) = chem
    chem.prune && !isnothing(tree) && !isempty(tree.living) && Trees.prune!(tree; lazy=true)
    nothing
end

# Brings the living cells up to date and sorts the tree, so that it can be analysed
update_tree!(::Nothing) = nothing
function update_tree!(tree::Tree)
    for cell in collect(keys(tree.living))
        Trees.cell_update!(tree, cell, get_curr_t(cell), get_data(cell))
    end

    sort!(tree)
end
