/*************************************************************
* @author   Triston Ruiseco
* @file     CLT.C
* @date     04/23/2021
* @brief    Script for generating a visualization of the
            manifestation of the central limit theorem
            as the number of observations averaged per
            sample is increased.
*************************************************************/

// C++ Standard Library Includes
#include <iostream>
#include <sstream>

// ROOT Includes
#include "TGraph.h"
#include "TVectorD.h"
#include "TFrame.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TColor.h"

// Global parameters and objects
const int numDims = 2;
const int numClusters = 4;
const int numPoints = 30;
const int maxIterations = 20;
double dimMax = 10.;
double dimMin = -10.;

// Directives and declarations for namespaces and namespace members
using namespace std;

// Definition for object representing data points in mixture
struct Point {

  // Storage variables
  double x, y;    // 2D coordinates
  int cluster;    // Nearest checked cluster
  double minDist; // Distance from nearest checked cluster

  // "Empty" constructor
  Point() : x(0.0), y(0.0), cluster(-1), minDist(__DBL_MAX__) {}

  // Constructor given data
  Point(double x, double y) : x(x), y(y), cluster(-1), minDist(__DBL_MAX__) {}

  // Distance from another Point
  double distance(Point p){
    return (p.x - x) * (p.x - x) + (p.y - y) * (p.y - y);
  }

};

// Begin primary program function
int main(int argc, char** argv){

  // Data generation parameter storage objects
  double means[numClusters][numDims];
  double stdevs[numClusters][numDims];
  vector<Point> data;

  // Randomly generate means and stdevs for each dimension of each cluster
  for(int c = 0; c < numClusters; c++){
    for(int d = 0; d < numDims; d++){
      means[c][d] = (dimMax - dimMin)*gRandom->Rndm() - dimMin;
      stdevs[c][d] = ((dimMax - dimMin)*gRandom->Rndm() - dimMin)/4.;
    }
  }

  // Generate and store data (only 2D, not yet generalized to numDims dimensions)
  for(int c = 0; c < numClusters; c++){
    for(int p = 0; p < numPoints; p++){
      data.push_back(Point(
        gRandom->Gaus(means[c][0], stdevs[c][0]),
        gRandom->Gaus(means[c][1], stdevs[c][1])));
    }
  }

  // Begin K-means clustering routine
  int dataSize = data.size();

  // Initialize centroids
  vector<Point> centroids;
  for(int i = 0; i < numClusters; i++){
    centroids.push_back(data[gRandom->Integer(dataSize)]);
  }

  // Refinement loop
  for(int l = 0; l < maxIterations; l++){

    // Assign points to a cluster
    for(int c = 0; c < numClusters; c++){
      for(int i = 0; i < dataSize; i++){
        double distance = centroids[c].distance(data[i]);
        if(distance < data[i].minDist){
          data[i].minDist = distance;
          data[i].cluster = c;
        }
      }
    }

    // Centroid refinement helper objects
    vector<int> nPoints;
    vector<double> sumX, sumY;

    // Initialize helper objects
    for(int i = 0; i < numClusters; i++){
      nPoints.push_back(0);
      sumX.push_back(0.);
      sumY.push_back(0.);
    }

    // Fill helper variables
    for(int i = 0; i < dataSize; i++){
      int c = data[i].cluster;
      nPoints[c] += 1;
      sumX[c] += data[i].x;
      sumY[c] += data[i].y;

      data[i].minDist = __DBL_MAX__;
    }

    // Calculate new centroids
    for(int c = 0; c < numClusters; c++){
      centroids[c].x = sumX[c] / double(nPoints[c]);
      centroids[c].y = sumY[c] / double(nPoints[c]);
    }

  }

  // Move results into visualization-friendly objects
  vector<double> xvals[numClusters];
  vector<double> yvals[numClusters];

  for(int i = 0; i < dataSize; i++){
    xvals[data[i].cluster].push_back(data[i].x);
    yvals[data[i].cluster].push_back(data[i].y);
  }

  // ROOT visualization objects and helper variables
  std::string namestr = "\0";
  TCanvas* canvas = new TCanvas("canvas1","canvas1", 700, 500);


  // ROOT visualization constant preferences
  canvas->SetGrid();
  canvas->DrawFrame(dimMin,dimMin,dimMax,dimMax);
  canvas->GetFrame()->SetBorderSize(12);

  // Array of groovy colors
  Int_t colors[9] = {TColor::GetColor("#f58582"),
                     TColor::GetColor("#dd9a44"),
                     TColor::GetColor("#7abd42"),
                     TColor::GetColor("#51c788"),
                     TColor::GetColor("#50c6ba"),
                     TColor::GetColor("#4ebaef"),
                     TColor::GetColor("#9b9afe"),
                     TColor::GetColor("#e876f0"),
                     TColor::GetColor("#fb74b7")};

  for(int c = 0; c < numClusters; c++){ // Loop across all clusters

    // Unique ROOT histogram name
    std::stringstream name;
    name << "cluster" << c;
    namestr = name.str();

    // ROOT graphs for visualization
    TGraph* graph = new TGraph(xvals[c].size(),&(xvals[c][0]),&(yvals[c][0]));
    graph->SetTitle(namestr.c_str());
    graph->SetName(namestr.c_str());

    // Visualization detail adjustments
    graph->SetMarkerColor(colors[2*c]);
    graph->SetMarkerStyle(20);

    // Cache finalized graph
    graph->Draw();
//    canvas->Update();
  }

  // Export final graphic
  canvas->SaveAs("KMeans.png");

  return 0;
}
