from ase.lattice.surface import fcc111_root

def check(root):
    try:
        slab = fcc111_root("Pt", root, (1,1,1), search_zone = (10,10))
        return True
    except:
        return False

# First test ability to generate existing confirguations in a range
valid = []
for root in range(10):
    if check(root):
        valid.append(root)

assert valid == [1, 3, 4, 7, 9]

# Now check if it fails when supplied bad arguments
failed = False
try:
    slab = fcc111_root("H", root, (1,1,1), search_zone = (10,10))
except:
    failed = True
assert failed 
