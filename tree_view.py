#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#  tree_view.py
#
#  Copyright 2014 Sébastien Gélis-Jeanvoine <sebastien@gelis.ch>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/>.
#

from ete2 import Tree, TreeStyle

tree = Tree("sequences.tree")
tree_style = TreeStyle()
tree_style.show_leaf_name = True
tree_style.show_branch_length = True
tree_style.show_branch_support = True
# tree_style.scale = 50
# tree_style.force_topology = True
# tree_style.mode = "c"
# tree_style.arc_start = 0
# tree_style.arc_span = 360
tree.render("sequences.png", tree_style=tree_style)
tree.show(tree_style=tree_style)