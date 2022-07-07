#include "Astar_searcher.h"
#include "node.h"
#include <Eigen/src/Core/Matrix.h>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <ros/assert.h>
#include <vector>


using namespace std;
using namespace Eigen;

void AstarPathFinder::initGridMap(double _resolution, Vector3d global_xyz_l,
                                  Vector3d global_xyz_u, int max_x_id,
                                  int max_y_id, int max_z_id) {
  gl_xl = global_xyz_l(0);
  gl_yl = global_xyz_l(1);
  gl_zl = global_xyz_l(2);

  gl_xu = global_xyz_u(0);
  gl_yu = global_xyz_u(1);
  gl_zu = global_xyz_u(2);

  GLX_SIZE = max_x_id;
  GLY_SIZE = max_y_id;
  GLZ_SIZE = max_z_id;
  GLYZ_SIZE = GLY_SIZE * GLZ_SIZE;
  GLXYZ_SIZE = GLX_SIZE * GLYZ_SIZE;

  resolution = _resolution;
  inv_resolution = 1.0 / _resolution;

  data = new uint8_t[GLXYZ_SIZE];
  memset(data, 0, GLXYZ_SIZE * sizeof(uint8_t));

  GridNodeMap = new GridNodePtr **[GLX_SIZE];
  for (int i = 0; i < GLX_SIZE; i++) {
    GridNodeMap[i] = new GridNodePtr *[GLY_SIZE];
    for (int j = 0; j < GLY_SIZE; j++) {
      GridNodeMap[i][j] = new GridNodePtr[GLZ_SIZE];
      for (int k = 0; k < GLZ_SIZE; k++) {
        Vector3i tmpIdx(i, j, k);
        Vector3d pos = gridIndex2coord(tmpIdx);
        GridNodeMap[i][j][k] = new GridNode(tmpIdx, pos);
      }
    }
  }
}

void AstarPathFinder::resetGrid(GridNodePtr ptr) {
  ptr->id = 0;
  ptr->cameFrom = NULL;
  ptr->gScore = inf;
  ptr->fScore = inf;
}

void AstarPathFinder::resetUsedGrids() {
  for (int i = 0; i < GLX_SIZE; i++)
    for (int j = 0; j < GLY_SIZE; j++)
      for (int k = 0; k < GLZ_SIZE; k++)
        resetGrid(GridNodeMap[i][j][k]);
}

void AstarPathFinder::setObs(const double coord_x, const double coord_y,
                             const double coord_z) {
  if (coord_x < gl_xl || coord_y < gl_yl || coord_z < gl_zl ||
      coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu)
    return;

  int idx_x = static_cast<int>((coord_x - gl_xl) * inv_resolution);
  int idx_y = static_cast<int>((coord_y - gl_yl) * inv_resolution);
  int idx_z = static_cast<int>((coord_z - gl_zl) * inv_resolution);

  if (idx_x == 0 || idx_y == 0 || idx_z == GLZ_SIZE || idx_x == GLX_SIZE ||
      idx_y == GLY_SIZE)
    data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
  else {
    data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
    data[(idx_x + 1) * GLYZ_SIZE + (idx_y + 1) * GLZ_SIZE + idx_z] = 1;
    data[(idx_x + 1) * GLYZ_SIZE + (idx_y - 1) * GLZ_SIZE + idx_z] = 1;
    data[(idx_x - 1) * GLYZ_SIZE + (idx_y + 1) * GLZ_SIZE + idx_z] = 1;
    data[(idx_x - 1) * GLYZ_SIZE + (idx_y - 1) * GLZ_SIZE + idx_z] = 1;
    data[(idx_x)*GLYZ_SIZE + (idx_y + 1) * GLZ_SIZE + idx_z] = 1;
    data[(idx_x)*GLYZ_SIZE + (idx_y - 1) * GLZ_SIZE + idx_z] = 1;
    data[(idx_x + 1) * GLYZ_SIZE + (idx_y)*GLZ_SIZE + idx_z] = 1;
    data[(idx_x - 1) * GLYZ_SIZE + (idx_y)*GLZ_SIZE + idx_z] = 1;
  }
}

vector<Vector3d> AstarPathFinder::getVisitedNodes() {
  vector<Vector3d> visited_nodes;
  for (int i = 0; i < GLX_SIZE; i++)
    for (int j = 0; j < GLY_SIZE; j++)
      for (int k = 0; k < GLZ_SIZE; k++) {
        // if(GridNodeMap[i][j][k]->id != 0) // visualize all nodes in open and
        // close list
        if (GridNodeMap[i][j][k]->id ==
            -1) // visualize nodes in close list only
          visited_nodes.push_back(GridNodeMap[i][j][k]->coord);
      }

  ROS_WARN("visited_nodes size : %d", visited_nodes.size());
  return visited_nodes;
}

Vector3d AstarPathFinder::gridIndex2coord(const Vector3i &index) {
  Vector3d pt;

  pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
  pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
  pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

  return pt;
}

Vector3i AstarPathFinder::coord2gridIndex(const Vector3d &pt) {
  Vector3i idx;
  idx << min(max(int((pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
      min(max(int((pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
      min(max(int((pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);

  return idx;
}

Eigen::Vector3d AstarPathFinder::coordRounding(const Eigen::Vector3d &coord) {
  return gridIndex2coord(coord2gridIndex(coord));
}

inline bool AstarPathFinder::isOccupied(const Eigen::Vector3i &index) const {
  return isOccupied(index(0), index(1), index(2));
}

inline bool AstarPathFinder::isFree(const Eigen::Vector3i &index) const {
  return isFree(index(0), index(1), index(2));
}

inline bool AstarPathFinder::isOccupied(const int &idx_x, const int &idx_y,
                                        const int &idx_z) const {
  return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE &&
          idx_z >= 0 && idx_z < GLZ_SIZE &&
          (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] == 1));
}

inline bool AstarPathFinder::isFree(const int &idx_x, const int &idx_y,
                                    const int &idx_z) const {
  return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE &&
          idx_z >= 0 && idx_z < GLZ_SIZE &&
          (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] < 1));
}

inline void AstarPathFinder::AstarGetSucc(GridNodePtr currentPtr,
                                          vector<GridNodePtr> &neighborPtrSets,
                                          vector<double> &edgeCostSets) {
  neighborPtrSets.clear();
  edgeCostSets.clear();
  Vector3i neighborIdx;
  expand_directions = {{-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 0, -1}, 
  {-1, 0, 0}, {-1, 0, 1}, {-1, 1, -1}, {-1, 1, 0}, {-1, 1, 1}, {0, -1, -1}, 
  {0, -1, 0}, {0, -1, 1}, {0, 0, -1}, {0, 0, 1}, {0, 1, -1}, {0, 1, 0}, 
  {0, 1, 1}, {1, -1, -1}, {1, -1, 0}, {1, -1, 1}, {1, 0, -1}, {1, 0, 0}, 
  {1, 0, 1}, {1, 1, -1}, {1, 1, 0}, {1, 1, 1}};
  expand_costs = {1.732, 1.414, 1.732, 1.414, 1, 1.414, 1.732, 1.414, 1.732,
									   1.414, 1, 1.414, 1, 1, 1.414, 1, 1.414,
									   1.732, 1.414, 1.732, 1.414, 1, 1.414, 1.732, 1.414, 1.732};
  expand_obstacles = {{-1,1,0},{1,1,0},{-1,-1,0},{1,-1,0}};
  
  for(int i = 0; i < expand_directions.size(); ++i)
  {
    neighborIdx = currentPtr->index + expand_directions[i];
    if (neighborIdx(0) < 0 || neighborIdx(0) >= GLX_SIZE ||
        neighborIdx(1) < 0 || neighborIdx(1) >= GLY_SIZE ||
        neighborIdx(2) < 0 || neighborIdx(2) >= GLZ_SIZE) {
      continue;
    }
    if (isOccupied(neighborIdx))
      continue;
    bool is_collison = false;
    for(auto &coner:expand_obstacles)
    {
      if(isOccupied(neighborIdx+coner)) 
      {
        is_collison = true;
        break;
      }
    }
    if (is_collison) continue;

    neighborPtrSets.push_back(
        GridNodeMap[neighborIdx(0)][neighborIdx(1)][neighborIdx(2)]);
    edgeCostSets.push_back(expand_costs[i]);
  }

 

  // for (int dx = -1; dx < 2; dx++) {
  //   for (int dy = -1; dy < 2; dy++) {
  //     for (int dz = -1; dz < 2; dz++) {

  //       if (dx == 0 && dy == 0 && dz == 0)
  //         continue;

  //       neighborIdx(0) = (currentPtr->index)(0) + dx;
  //       neighborIdx(1) = (currentPtr->index)(1) + dy;
  //       neighborIdx(2) = (currentPtr->index)(2) + dz;

  //       if (neighborIdx(0) < 0 || neighborIdx(0) >= GLX_SIZE ||
  //           neighborIdx(1) < 0 || neighborIdx(1) >= GLY_SIZE ||
  //           neighborIdx(2) < 0 || neighborIdx(2) >= GLZ_SIZE) {
  //         continue;
  //       }
  //       if (isOccupied(dx,dy,dz)) continue;

  //       neighborPtrSets.push_back(
  //           GridNodeMap[neighborIdx(0)][neighborIdx(1)][neighborIdx(2)]);
  //       edgeCostSets.push_back(sqrt(dx * dx + dy * dy + dz * dz));
  //     }
  //   }
  // }
}

double AstarPathFinder::getHeu(GridNodePtr node1, GridNodePtr node2) {
  // using digonal distance and one type of tie_breaker.
  double h = inf;
  double dx = fabs(node1->coord.x() - node2->coord.x());
  double dy = fabs(node1->coord.y() - node2->coord.y());
  double dz = fabs(node1->coord.z() - node2->coord.z());

  h = sqrt(dx * dx + dy * dy + dz * dz);

  double kTieBreaker = 0.1;

  return h * (1+kTieBreaker);
}

void AstarPathFinder::AstarGraphSearch(Vector3d start_pt, Vector3d end_pt) {
  ros::Time time_1 = ros::Time::now();

  // index of start_point and end_point
  Vector3i start_idx = coord2gridIndex(start_pt);
  Vector3i end_idx = coord2gridIndex(end_pt);
  goalIdx = end_idx;

  // position of start_point and end_point
  start_pt = gridIndex2coord(start_idx);
  end_pt = gridIndex2coord(end_idx);

  // Initialize the pointers of struct GridNode which represent start node and
  // goal node
  GridNodePtr startPtr = new GridNode(start_idx, start_pt);
  GridNodePtr endPtr = new GridNode(end_idx, end_pt);

  // openSet is the open_list implemented through multimap in STL library
  openSet.clear();
  // currentPtr represents the node with lowest f(n) in the open_list
  GridNodePtr currentPtr = NULL;
  GridNodePtr neighborPtr = NULL;

  // put start node in open set
  startPtr->gScore = 0;
  /**
   *
   * STEP 1.1:  finish the AstarPathFinder::getHeu
   *
   * **/
  startPtr->fScore = getHeu(startPtr, endPtr);  // 使用了MANHANTT距离

  startPtr->id = 1;
  startPtr->coord = start_pt;
  startPtr->cameFrom = nullptr;
  openSet.insert(make_pair(startPtr->fScore, startPtr));

  /**
   *
   * STEP 1.2:  some else preparatory works which should be done before while
   * loop
   *
   * **/

  // double tentative_gScore;
  vector<GridNodePtr> neighborPtrSets;
  vector<double> edgeCostSets;

  /**
   *
   * STEP 1.3:  finish the loop
   *
   * **/
  while (!openSet.empty()) {


    auto iter = openSet.begin(); 
        currentPtr = (*iter).second;
        currentPtr->id = -1;
        openSet.erase((*iter).first);

        // if the current node is the goal 
        if( currentPtr->index == goalIdx ){
            ros::Time time_2 = ros::Time::now();
            terminatePtr = currentPtr;
            ROS_INFO("[A*]{sucess}  Time in A*  is %f ms, path cost if %f m", (time_2 - time_1).toSec() * 1000.0, currentPtr->gScore * resolution );            
            return;
        }
        //get the succetion
        AstarGetSucc(currentPtr, neighborPtrSets, edgeCostSets);
        for(int i = 0; i < (int)neighborPtrSets.size(); i++){
            /*
            *
            *
            Judge if the neigbors have been expanded
            please write your code below
            
            IMPORTANT NOTE!!!
            neighborPtrSets[i]->id = -1 : expanded, equal to this node is in close set
            neighborPtrSets[i]->id = 1 : unexpanded, equal to this node is in open set
            *        
            */
            neighborPtr = neighborPtrSets.at(i);
            double gScore = currentPtr->gScore + edgeCostSets.at(i);
            double hScore = getHeu(neighborPtr, endPtr);
            double fScore = gScore + hScore;

            if(neighborPtr -> id == 0){ //discover a new node, which is not in the closed set and open set

                neighborPtr->gScore = gScore;
                neighborPtr->fScore = fScore;
                neighborPtr->id = 1; // OPEN set
                neighborPtr->cameFrom = currentPtr;
                openSet.insert(std::make_pair(neighborPtr->fScore,neighborPtr));
                
                continue;
            }
            else if(neighborPtr->id == 1){ //this node is in open set and need to judge if it needs to update, the "0" should be deleted when you are coding

                if(neighborPtr->gScore > gScore) {
                    
                    neighborPtr->gScore = gScore;
                    neighborPtr->fScore = fScore;
                    neighborPtr->cameFrom = currentPtr;
                    openSet.insert(std::make_pair(neighborPtr->fScore,neighborPtr)); 
                }

                continue;
            }
            else{//this node is in closed set

                continue;
            }
        }



  }

  // if search fails
  ros::Time time_2 = ros::Time::now();
  if ((time_2 - time_1).toSec() > 0.1)
    ROS_WARN("Time consume in Astar path finding is %f",
             (time_2 - time_1).toSec());
}

vector<Vector3d> AstarPathFinder::getPath() {
  vector<Vector3d> path;
  // vector<GridNodePtr> gridPath;

  /**
   *
   * STEP 1.4:  trace back the found path
   *
   * **/
    path.clear();

    while (terminatePtr)
    {
        path.push_back(terminatePtr->coord);
        terminatePtr = terminatePtr->cameFrom;
    }
           
    reverse(path.begin(),path.end());

  return path;
}



vector<Vector3d> AstarPathFinder::pathSimplify(const vector<Vector3d> &path,
                                               double path_resolution) {
  // implement the RDP algorithm
  ROS_INFO("in function pathSimplify - init");

  vector<Vector3d> subPath;
  subPath.clear();
  if(path.size() < 2) return subPath;
  subPath.push_back(path.front());

  auto start_iter = path.begin();
  auto inter_iter = path.end() - 1;
  // std::cout << "start pt before while loop:" << *start_iter << std::endl;
  // std::cout << "inter pt before while loop:" << *inter_iter << std::endl;
    
  double max_distance_threshold = path_resolution * 2.0; // !搜索步长相关
  while (start_iter != inter_iter) {
    Vector3d vec_b = *start_iter - *inter_iter;
    
    double segment_length = vec_b.norm();
    segment_length =
        segment_length < 1e-6 ? 0.00001 : segment_length; // 避免被除数为0
    
    // ROS_INFO("segment_length : %f",segment_length);
    bool isfound = false;
    for (auto it = start_iter; it != inter_iter; it++) {
      Vector3d vec_a = *it - *start_iter;
      // if (vec_a.norm() < 0.001) continue;
      double area = vec_b.cross(vec_a).norm();
      double vertical_dist = area / segment_length;
      // ROS_INFO("area : %f , vertical_dist : %f ",segment_length,vertical_dist);
      if (vertical_dist > max_distance_threshold) {
        inter_iter = it;
        isfound = true;
        break;
        }
    }
    // std::cout << "start pt :" << *start_iter << std::endl;
    // std::cout << "inter pt :" << *inter_iter << std::endl;
    if (!isfound) {
      subPath.push_back(*inter_iter);
      start_iter = inter_iter;
      inter_iter = path.end()-1;
    }
  }
  ROS_INFO("subPath size : %d", subPath.size());
  return subPath;
}

Vector3d AstarPathFinder::getPosPoly(MatrixXd polyCoeff, int k, double t) {
  Vector3d ret;
  int _poly_num1D = (int)polyCoeff.cols() / 3;
  for (int dim = 0; dim < 3; dim++) {
    VectorXd coeff = (polyCoeff.row(k)).segment(dim * _poly_num1D, _poly_num1D);
    VectorXd time = VectorXd::Zero(_poly_num1D);

    for (int j = 0; j < _poly_num1D; j++)
      if (j == 0)
        time(j) = 1.0;
      else
        time(j) = pow(t, j);

    ret(dim) = coeff.dot(time);
    // cout << "dim:" << dim << " coeff:" << coeff << endl;
  }

  return ret;
}

int AstarPathFinder::safeCheck(MatrixXd polyCoeff, VectorXd time) {
  int unsafe_segment = -1; //-1 -> the whole trajectory is safe
  /**
   *
   * STEP 3.3:  finish the safeCheck()
   *
   * **/
  for (int k = 0; k < time.size(); ++k) {
    for (double t = 0; t <= time[k]; t += 0.1) {
      auto pos = getPosPoly(polyCoeff, k, t);  // 无人机不能只当作一个质点
      if (isOccupied(coord2gridIndex(pos))) {
        unsafe_segment = k;
        break;
      }
    }
  }
  

  return unsafe_segment;
}