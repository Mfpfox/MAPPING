from pymol import cmd, CmdException


# def alphaToAll(sel, col="b",forceRebuild=True):

def alphaToAll(sel, col="b",forceRebuild=True):
	"""
	alphaToAll -- expand any property of the alpha carbons to all atoms in the residue
	
	PARAMS
		sel
			The selection or object (include "*" for all objects) to operate on.  This will
			read the alpha carbon's "COL" property and expand that down to all atoms in the
			given residues in sel.
			
		col
			Any valid PyMOL property.  For example, 'b' or 'q' or 'color', etc.
			DEFAULT: b, for B-factor as we most often overwrite that column
			
		forceRebuild
			If a color changes, this will rebuild the objects for us.  For speed, we
			DEFAULT this to False.
			
	RETURNS
		None, just epxands the properties as dicsussed above.
		
	NOTES
		This script is useful for instances where we update the b-factor column for the
		alpha carbons and then want to expand that color to all atoms for smoother looking
		images.
		
	AUTHOR:
		Jason Vertrees, 2009.
			
	"""
	col = '(' + ','.join(col.split()) + ')'
	space = {'props': dict()}
	cmd.iterate('byca (%s)' % (sel), 'props[model,segi,chain,resi] = ' + col, space=space)
	cmd.alter(sel, col + ' = props.get((model,segi,chain,resi), ' + col + ')', space=space)
	if forceRebuild != False:
		cmd.rebuild(sel)

cmd.extend("alphaToAll", alphaToAll)