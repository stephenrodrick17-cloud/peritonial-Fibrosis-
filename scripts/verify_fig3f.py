import re

svg = open("results/revision/matrisome_string_ppi/Fig3F_STRING_44_cleaned.svg", "r").read()

# Count node containers
nodes = re.findall(r"class='nwnodecontainer'", svg)
print("Node containers:", len(nodes))

# Extract gene labels
labels = re.findall(r"class='nwnodelabel'>([A-Za-z0-9]+)<", svg)
print("Labels found:", len(labels))
print("All labels:", sorted(labels))
print()
print("NAT6 in SVG:", "NAT6" in svg)
print("MBP in SVG:", "MBP" in svg)
print("CCL5 in SVG:", "CCL5" in svg)

# Count edges
edges = re.findall(r"class='nwlinkwrapper'", svg)
print("Edge groups:", len(edges))

# Check degree of CCL5 by counting edges that connect to it
# Find CCL5 node id first
ccl5_refs = re.findall(r"id='(edge\.[^']*)'", svg)
print("Total edge ids:", len(ccl5_refs))
