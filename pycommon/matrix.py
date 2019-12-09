###### bin N-column matrix to 2 column with counts
### using python
def matrix_contract(mat, cx, cy, cz):
    #c_x = mat[:][cx-1]
    #c_y = mat[:][cy-1]
    #c_z = mat[:][cz-1]
    x = []
    y = []
    z = []
    xyz = []
    h, w = mat.shape
    for i in range(h):
        if mat[i,cx-1] in x and mat[i,cy-1] in y:
            index = 
            
        
