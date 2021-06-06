#include <iostream>
#include <vector>
#include <cstdio>
#include <cstring>

struct coordinate {
  double x;
  double y;
  double z;
};
double distance_squared(const coordinate & r1, const coordinate & r2) {
  double sum, delta;

  sum = 0;
  delta = r1.x - r2.x;  sum += delta * delta;
  delta = r1.y - r2.y;  sum += delta * delta;
  delta = r1.z - r2.z;  sum += delta * delta;
  return sum;

}
void  get_surface(const char * protein_name, const char * output_name, int num_surface_points) {
  
  int     num_neighbors_buried = 10;
  double  distance_neighbor_2  = 4*4;   // squared distance between atom and its neighbor

  // Loading protein
  FILE * f = fopen(protein_name, "r");
  if (f == NULL) {  
    std::cout << "Unable to open file " << protein_name << "\n";
  }
  std::vector<coordinate>     coords;
  std::vector<std::string>    pdb_record;
  char   atom_name[5];
  char   residue_name[4];
  char   buff [256];
  while (fgets(buff, sizeof buff, f) != NULL) {
    if ((strstr(buff, "ATOM")!= NULL || strstr(buff, "HETATM")!= NULL) && strlen(buff) > 57) {
      
      // Hydrogens are discarded:
      sscanf(buff+11, "%4s", atom_name);
      if (atom_name[0] == 'H' || (atom_name[0] > 0 && atom_name[0] <= '9' && atom_name[1] == 'H')) continue;

      // Water molecules are discarded:
      sscanf(buff+17, "%3s", residue_name);
      if (strcmp(residue_name, "HOH")==0 || strcmp(residue_name, "WAT")==0) continue;

      coordinate r;
      //double r[3];
      sscanf(buff+31, "%lf %lf %lf", &r.x, &r.y, &r.z);
      coords.push_back(r);
      pdb_record.push_back(buff);

    }
  }
  fclose(f);

  // Calculating, how many neighbors each atom has
  std::vector<int> num_neighbors(coords.size());
  for (int i=0; i<coords.size(); i++) {
    num_neighbors[i] = 0;
    for (int j=0; j<coords.size(); j++) {
      double r2 = distance_squared(coords[i], coords[j]);
      if (r2 < distance_neighbor_2)  
        num_neighbors[i]++;
    }
  }

  // Getting rid of buried atoms
  std::vector<coordinate>   tmp_coords;
  std::vector<std::string>  tmp_pdb_record;
  for (int i=0; i<coords.size(); i++) {
    if (num_neighbors[i] < num_neighbors_buried) {
      tmp_coords.push_back(coords[i]);
      tmp_pdb_record.push_back(pdb_record[i]);
    }
  }
  coords = tmp_coords; 
  pdb_record = tmp_pdb_record;
  
  // Selecting atoms for output
  int natoms = coords.size();
  if (num_surface_points > coords.size()) num_surface_points = coords.size();
  std::vector<double> min_distance(natoms); // Distance between i'th atom and nearest taken atom
  std::vector<bool> atom_taken(natoms);     // Was i'th atom taken 
  std::vector<int>    surface_atoms;        // List of taken atoms
  surface_atoms.push_back(0);
  for (int i=0; i<natoms; i++) {
    min_distance[i] = distance_squared(coords[0], coords[i]);
    atom_taken[i] = false;
  }
  atom_taken[0] = true;
  for (int n=0; n<num_surface_points; n++) {
    double maxmin_dist = 0;
    int    best_atom   = -1;
    for (int i=0; i<natoms; i++) if (!atom_taken[i] && min_distance[i] > maxmin_dist)  {
      maxmin_dist = min_distance[i];
      best_atom   = i;
    }
    if (best_atom < 0) continue;
    surface_atoms.push_back(best_atom);
    atom_taken[best_atom] = true;

    for (int i=0; i<natoms; i++) if (!atom_taken[i]) {
      double new_distance = distance_squared(coords[i], coords[best_atom]);
      if (new_distance < min_distance[i]) {
        min_distance[i] = new_distance;
      }
    }
  }

  // Saving output in pdb format  
  f = fopen(output_name,"w");
  for (int i=0; i<natoms; i++) if (atom_taken[i]) {
    fprintf(f,"%s", pdb_record[i].c_str());
  }
  fclose(f);

}

int   main(int argc, char* argv[]) {	
  std::cout << "usage: surface input.pdb output.pdb num\n";
  if (argc < 4) return 1;
  int num_points; 
  sscanf(argv[3],"%i",&num_points);
  get_surface(argv[1], argv[2], num_points );
  
}