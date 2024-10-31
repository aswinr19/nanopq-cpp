#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

void Vq(std::vector<std::vector<std::vector<double>>> &codewords,
        std::vector<std::vector<double>> &subVectors,
        std::vector<std::vector<int>> &codes, int dim1, int dim2, int dim3,
        int dim4, int dim5, int ds, int fixed) {

  std::vector<double> shortestDistIdxs(dim4, 0); // 2000 x 1
  double currentDist, shortestDist;
  int shortestIdx = 0;

  for (int i = 0; i < dim4; i++) {
    shortestDist = 1000;
    currentDist = 1000;

    for (int j = 0; j < dim2; j++) {

      double sum = 0;
      for (int k = 0; k < dim3; k++) {
        sum += pow((codewords[fixed][j][k] - subVectors[i][k]), 2);
      }
      currentDist = sqrt(sum);
      if (currentDist <= shortestDist) {
        shortestDist = currentDist;
        shortestIdx = j;
      }
    }
    codes[i][fixed] = shortestIdx;
  }
}

std::vector<std::vector<int>>
Encode(std::vector<std::vector<std::vector<double>>> &codewords,
       std::vector<std::vector<double>> &vectors, int dim1, int dim2, int dim3,
       int dim4, int dim5, int m) {

  int ds = dim5 / m;
  int p, q;

  std::vector<std::vector<int>> codes(dim4, std::vector<int>(m, 0)); // 2000 x 8

  for (int i = 0; i < m; i++) {

    std::vector<std::vector<double>> subVecs(dim4, std::vector<double>(ds, 0));

    for (int j = 0; j < dim4; j++) {
      for (int k = 0; k < ds; k++) {
        for (int l = i * ds; l < (i + 1) * ds; l++) {
          subVecs[j][k] = vectors[j][l];
        }
      }
    }
    Vq(codewords, subVecs, codes, dim1, dim2, dim3, dim4, dim5, ds, i);
  }
  return codes;
}

std::vector<std::vector<double>> ReadVectors(std::string path, int dim1,
                                             int dim2) {

  std::vector<std::vector<double>> vectors(dim1, std::vector<double>(dim2, 0));
  std::ifstream vectorFile(path);

  if (!vectorFile.is_open())
    throw std::runtime_error("Could not open file!\n");

  std::vector<std::vector<double>> vecs(dim1, std::vector<double>(dim2, 0));
  std::string line;
  int count = 0;

  if (vectorFile.good()) {
    int firstIdx = 0;
    int secondIdx = 0;

    while (std::getline(vectorFile, line) && firstIdx < dim1) {
      std::stringstream ss(line);
      std::string cell;

      secondIdx = 0;
      while (std::getline(ss, cell, ',')) {
        count++;
        vectors[firstIdx][secondIdx] = stod(cell);
        secondIdx++;
      }
      firstIdx++;
    }
  }

  return vectors;
}

std::vector<std::vector<std::vector<double>>>
ReadCsv(std::string path, int dim1, int dim2, int dim3) {

  std::ifstream csvFile(path);
  std::vector<std::vector<std::vector<double>>> codewords(
      dim1,
      std::vector<std::vector<double>>(dim2, std::vector<double>(dim3, 0)));

  if (!csvFile.is_open())
    throw std::runtime_error("Could not open file!\n");

  std::string line;

  if (csvFile.good()) {
    int firstIdx = 0;
    int secondIdx = 0;
    int thirdIdx = 0;

    while (std::getline(csvFile, line) && firstIdx < dim1) {
      std::stringstream ss(line);
      std::string cell;

      secondIdx = 0;
      while (std::getline(ss, cell, ',') &&
             secondIdx * thirdIdx < dim2 * dim3) {

        codewords[firstIdx][secondIdx][thirdIdx] = stod(cell);
        if (thirdIdx < dim3 - 1) {
          thirdIdx++;
        } else {
          secondIdx++;
          thirdIdx = 0;
        }
      }
      firstIdx++;
    }
  }
  return codewords;
}

std::vector<int> ReadDimensions(std::string path) {
  std::ifstream shapeFile(path);

  std::vector<int> dimensions(10, 0);

  if (!shapeFile.is_open())
    throw std::runtime_error("Could not open shape file!\n");

  int count = 0;

  std::string line;
  if (shapeFile.good()) {
    while (std::getline(shapeFile, line)) {
      double value;
      std::istringstream iss(line);
      // convert exponential notation into decimal values
      if ((iss >> value)) {
        dimensions[count] = value;
        count++;
      }
    }
  }

  return dimensions;
}

int main(int argc, char **argv) {

  std::vector<int> dim(10, 0);

  dim = ReadDimensions(argv[1]);

  std::vector<std::vector<int>> codes(2000, std::vector<int>(8, 0));
  std::vector<std::vector<double>> vecs(dim[3], std::vector<double>(dim[4], 0));

  std::vector<std::vector<std::vector<double>>> codewords(
      dim[0],
      std::vector<std::vector<double>>(dim[1], std::vector<double>(dim[2], 0)));

  codewords = ReadCsv(argv[2], dim[0], dim[1], dim[2]);
  vecs = ReadVectors(argv[3], dim[3], dim[4]);
  codes =
      Encode(codewords, vecs, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]);

  for (int i = 0; i < 2000; i++) {
    for (int j = 0; j < 8; j++) {
      std::cout << codes[i][j] << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}
