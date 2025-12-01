#include <iterator>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include <iostream>
#include <tuple>
#include <list>
#include <cstring>
#include "local_search.hpp"

namespace py = pybind11;

template <typename T> 
T tourDistance(const std::vector<int>& tour, const T* disMat, int nnode) {
    T totalDistance = 0.0;
    for (size_t i = 0; i < tour.size(); ++i) {
        totalDistance += disMat[tour[i] * nnode +  tour[(i + 1) % nnode]];
    }
    return totalDistance;
}
// 2-opt swap function
void twoOptSwap(std::vector<int>& tour, int i, int k) {
    reverse(tour.begin() + i + 1, tour.begin() + k + 1);
}


// Function to find the nearest neighbor and update the tour
int findNearestNeighbor(int current_city, std::vector<bool>& visited, double* dist){
    int nearestCity = -1;
    double minDistance = 1e200;

    for(int i =0; i < visited.size(); i++){
        if (!visited[i] && i != current_city) {
            if(dist[i] < minDistance){
                minDistance = dist[i];
                nearestCity = i;
            }
        }
    }
    return nearestCity;
}


// Function to solve the TSP using the Nearest Neighbor algorithm

std::tuple<bool, py::array_t<int>>  decodeNNTSP(py::array_t<double> distMat, py::array_t<int> assignment, int start_node) {

    int nnode = assignment.size();
    std::vector<bool> visited(nnode, false);
    auto res = py::array_t<int>(nnode);
    
    bool assign_valid = true;



    py::buffer_info assignment_buf = assignment.request();
    py::buffer_info distMat_buf = distMat.request();
    py::buffer_info res_buf = res.request();

    double *ptr_distMat = static_cast<double *>(distMat_buf.ptr);
    int *ptr_assignment  = static_cast<int *>(assignment_buf.ptr);
    int *ptr_res = new int[nnode];
    int *ptr_res_ = static_cast<int *>(res_buf.ptr);

    int currentCity = start_node;
    int tour_len = 0;
    ptr_res[tour_len] = start_node;
    tour_len++;
    visited[start_node] = true;
    for(int i = 1; i < nnode; ++i){
        int next_city = -1;
        if(visited[ptr_assignment[currentCity]] == false){
            next_city = ptr_assignment[currentCity];
        }
        else{
            assign_valid = false;
            next_city = findNearestNeighbor(currentCity, visited, ptr_distMat + nnode * currentCity);
        }
        ptr_res[i] = next_city;
        currentCity = next_city;
        visited[currentCity] = true;
    }

    for(int i_ = 0; i_ < nnode; i_++){
        int i = ptr_res[i_];
        int j = ptr_res[(i_ + 1) % nnode];
        ptr_res_[i] = j;
    }

    delete []ptr_res;
    return std::make_tuple(assign_valid, res);
}


// 2-opt algorithm
template<typename T> 
double twoOpt(T* distMat, int nnode, std::vector<int>& tour) {
    // Start with an initial tour (e.g., sequential order)
    
    double bestDistance = tourDistance(tour, distMat, nnode);
    bool improved = true;

    

    while (improved) {
        improved = false;
        for (size_t i = 0; i < nnode - 1; ++i) {
            for (size_t k = i + 2; k < nnode; ++k) {
                std::vector<int> newTour = tour;
                twoOptSwap(newTour, i, k);
                double newDistance = tourDistance(newTour, distMat, nnode);

                if (newDistance < bestDistance) {
                    tour = newTour;
                    bestDistance = newDistance;
                    improved = true;
                    break;
                }
            }
            if(improved){
                break;
            }
        }
    }

    return bestDistance;
}


py::array_t<int> local_search_2opt(py::array_t<double> C, py::array_t<int> Path, bool sym=false){
    auto res = py::array_t<int>(Path.size());
    int d = Path.size();

    std::vector<int> tour;

    py::buffer_info Path_buf = Path.request();
    py::buffer_info res_buf = res.request();
    py::buffer_info C_buf = C.request();

    int *ptr_res = static_cast<int *>(res_buf.ptr);
    int *ptr_Path = static_cast<int *>(Path_buf.ptr);
    double *ptr_C = static_cast<double *>(C_buf.ptr);
    
    int cnode = 0;

    for(int i = 0; i < d; i++){
        tour.push_back(cnode);
        cnode = ptr_Path[cnode];
    }

    while(twoOptMove<double>(tour, ptr_C, Path.size(), sym));
    //twoOpt<double>(ptr_C, Path.size(), tour);
    for(int i = 0; i < d; i++){
        ptr_res[tour[i]] = tour[(i + 1) % Path.size()];
    }
    
    return res;
}




py::array_t<double> partial_tsp(py::array_t<double> C, py::array_t<int> cols, double max_v=1e20){
    auto C_hat = py::array_t<double>(C.size());
    std::vector<int> visited(cols.size());
    for(int i = 0; i < cols.size(); i++) visited[i] = 0;
    int vcnt = 0;
    
    py::buffer_info cols_buf = cols.request();
    py::buffer_info C_hat_buf = C_hat.request();
    py::buffer_info C_buf = C.request();
    
    int *ptr_cols = static_cast<int *>(cols_buf.ptr);
    double *ptr_C_hat = static_cast<double *>(C_hat_buf.ptr);
    double *ptr_C = static_cast<double *>(C_buf.ptr);
    for(int i = 0; i < C.size(); i++){
	ptr_C_hat[i] = max_v;
    }
    
    std::vector<std::vector<int>> tours;
    while(1){
	std::vector<int> tour;
	for(int j = 0; j < visited.size(); j++){
	    if(visited[j] == 0) {
		tour.push_back(j);
		visited[j] = 1;
		vcnt += 1;
		break;
	    }
	}
	if(vcnt == cols.size()) break;
	while(1){
	    int j = ptr_cols[tour[tour.size() - 1]];
	    if(visited[j]) break;
	    visited[j] = 1;
	    vcnt ++;
	    tour.push_back(j);
	}
	tours.push_back(tour);
	for(int i = 0; i < tour.size(); i++){
	    int i_ = tour[i];
	    int j_ = tour[(i + 1) % tour.size()];
	    ptr_C_hat[i_ * cols.size() + j_] = ptr_C[i_ * cols.size() + j_];
	}
	if(vcnt == cols.size()) break;
    }
    for(int i = 0; i < tours.size(); i++){
	for(int j = i + 1; j < tours.size(); j++){
	    for(auto i_: tours[i]){
		for(auto j_: tours[j]){
		    ptr_C_hat[i_ * cols.size() + j_] = ptr_C[i_ * cols.size() + j_];
		    ptr_C_hat[j_ * cols.size() + i_] = ptr_C[j_ * cols.size() + i_];
		}
	    }
	}
    }
    return C_hat;
}

py::array_t<int> inv_perm(py::array_t<int>cols){
    if (cols.ndim() != 1)
        throw std::runtime_error("Number of dimensions must be one");
    auto rows = py::array_t<int>(cols.size());


    py::buffer_info rows_buf = rows.request();
    py::buffer_info cols_buf = cols.request();

    int *ptr_cols = static_cast<int*>(cols_buf.ptr);
    int *ptr_rows = static_cast<int*>(rows_buf.ptr);

    for(int i = 0; i < rows.size(); i++){
        ptr_rows[ptr_cols[i]] = i;
    }
    return rows;
}

std::tuple<double, py::array_t<double>> partial_dag_cons(py::array_t<int> rows, py::array_t<int> cols, double scale){
    if (rows.ndim() != 1 || cols.ndim() != 1)
        throw std::runtime_error("Number of dimensions must be one");

    if (rows.size() != cols.size())
        throw std::runtime_error("Input shapes must match");

    auto result = py::array_t<double>(rows.size() * rows.size());
    double hvalue=0.0; 
    py::buffer_info result_buf = result.request();
    py::buffer_info rows_buf = rows.request();
    py::buffer_info cols_buf = cols.request();
    
    int *ptr_row = static_cast<int *>(rows_buf.ptr);
    int *ptr_col = static_cast<int *>(cols_buf.ptr);
    double *ptr_grad = static_cast<double *>(result_buf.ptr);

    int *ptr_inv_row = new int[rows.size()];
    int *ptr_power_perm = new int[rows.size()];
    int *ptr_power_perm_buf = new int[rows.size()];
    
    for(int i = 0; i < rows.size(); i++){
	ptr_inv_row[ptr_col[i]] = i;

	ptr_power_perm[i] = ptr_col[i];
	for(int j = 0; j < rows.size(); j++){
	    ptr_grad[i * rows.size() + j] = 0;
	}
	if(ptr_power_perm[i] == i) hvalue += scale;
    }
    double cscale = scale;
    for(int i = 1; i < rows.size() - 1; i++)
    {
	for(int j = 0; j < rows.size(); j++){
	    ptr_grad[ptr_power_perm[j] * rows.size() + j] += (i + 1) * cscale * scale;
#ifdef DEBUG
	    if(i==1)
		std::cout << ptr_inv_row[j] << " ";
#endif
	}
#ifdef DEBUG
	if(i==1) std::cout << std::endl;
#endif 

	cscale *= scale;
	
	for(int j=0; j < rows.size(); j++){
	    ptr_power_perm_buf[j] = ptr_col[ptr_power_perm[j]];
	}
	for(int j=0; j < rows.size(); j++){
	    ptr_power_perm[j] = ptr_power_perm_buf[j];
#ifdef DEBUG
	    std::cout << ptr_power_perm[j] << " ";
#endif
	    if(ptr_power_perm[j] == j) hvalue += scale;

	}
#ifdef DEBUG
	std::cout << std::endl;
#endif 
    }


    delete []ptr_inv_row;
    delete []ptr_power_perm;
    delete []ptr_power_perm_buf;

    return std::make_tuple(hvalue, result);
}

PYBIND11_MODULE(partial_dag_model, m) {
    m.def("partial_dag_cons", &partial_dag_cons, "return the partial dag constraints jac for a TSP problem");
    m.def("partial_tsp", &partial_tsp, "return the partial TSP problem");
    m.def("local_search_2opt", &local_search_2opt, "local search use 2opt");
    m.def("decodeNNTSP", &decodeNNTSP, "decode from assignment using NN");
}
