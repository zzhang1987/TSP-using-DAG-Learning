import numpy as np

def decompose_into_powers_of_two(n):
    """
    Decomposes an integer into the sum of powers of 2.
    
    Args:
        n (int): The integer to decompose. Must be >= 0.

    Returns:
        List[int]: A list of powers of 2 that sum to `n`.
    """
    if n < 0:
        raise ValueError("Only non-negative integers are supported.")

    powers = []
    power = 0
    while n > 0:
        if n % 2 == 1:
            powers.append(power)
        n = n // 2
        power += 1

    return powers


def load_single_problem_from_file(filename, node_cnt=None, scaler=1e6):

    ################################
    # "tmat" type
    ################################
    if node_cnt is not None:
        problem = torch.empty(size=(node_cnt, node_cnt), dtype=torch.long)
    # shape: (node, node)

    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except Exception as err:
        print(str(err))

    all_mat = []
    line_cnt = 0
    for line in lines:
        if line == '\n':
            continue
        linedata = line.split()

        if linedata[0].startswith(('NAME', 'COMMENT', 'TYPE', 'DIMENSION', 'EDGE_WEIGHT_TYPE', 'EDGE_WEIGHT_FORMAT', 'EDGE_WEIGHT_SECTION', 'EOF')):
            continue

        integer_map = map(int, linedata)
        integer_list = list(integer_map)

        #if node_cnt is None:
        #    node_cnt = len(integer_list)
        #    problem = torch.empty(size=(node_cnt, node_cnt), dtype=torch.long)

        #problem[line_cnt] = torch.tensor(integer_list, dtype=torch.long)
        all_mat += integer_list
        line_cnt += 1

    node_cnt = int(np.sqrt(len(all_mat)))
    problem = torch.from_numpy(np.asarray(all_mat).reshape(node_cnt, node_cnt))
                   
    # Diagonals to 0
    problem[torch.arange(node_cnt), torch.arange(node_cnt)] = 0

    # Scale
    scaled_problem = problem.float() / scaler

    return scaled_problem.double().numpy()

    

if __name__ == '__main__':
    import numpy as np 
    import logging
    import sys

    root = logging.getLogger()
    root.setLevel(logging.DEBUG)

    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    root.addHandler(handler)
    def testity(num, dec):
        logging.debug("test {} and decomposition {}".format(num, dec))
        assert(np.sum([2 ** i for i in dec]).item() == num)

        return True
            
    
    for i in range(10):
        j = np.random.randint(0, 65535)
        dec = decompose_into_powers_of_two(j)
        print(testity(j, dec))
