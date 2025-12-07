import os 
import numpy as np

def create_lkh_file(C: np.ndarray, filename: str, scale: int=1000000):
    nnode = C.shape[0]
    with open(filename, 'w') as f:
        f.write("NAME : temp\n")
        f.write("TYPE : ATSP\n")
        f.write(f"DIMENSION : {nnode}\n")
        f.write("EDGE_WEIGHT_TYPE : EXPLICIT\n")
        f.write("EDGE_WEIGHT_FORMAT : FULL_MATRIX\n")
        f.write("EDGE_WEIGHT_SECTION\n")
        for i_ in range(C.shape[0]):
            all_str=""
            for j in range(C.shape[0]):
                number = int(C[i_, j] * scale)
                desired_length = 11
                if j!= C.shape[0] - 1:
                    cstr = f"{number:<{desired_length}} "
                else:
                    cstr = f"{number:<{desired_length}}"
                all_str = all_str + cstr
            f.writelines(all_str + "\n")
        f.write("EOF\n")