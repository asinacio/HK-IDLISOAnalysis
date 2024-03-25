//*********************************************
// Store the locations of UKLI injectors inside
// Super-Kamiokande, for easy access in filling
// event trees.
//
// s.j.jenkins@liverpool.ac.uk
//********************************************

#ifndef SKINJECTORLOCATIONS_H
#define SKINJECTORLOCATIONS_H

#include <vector>
#include <string>

std::vector<double> GetInjectorLocation(std::string inj){

  std::vector<double> vtx(3);

  //All x and y coordinates are common
  vtx[0] = 1490.73;
  vtx[1] = 768.14;

  if(inj == "B1")      vtx[2] = 1223.9;
  else if(inj == "B2") vtx[2] = 742.3;
  else if(inj == "B3") vtx[2] = -200.3;
  else if(inj == "B4") vtx[2] = -747.1;
  else if(inj == "B5") vtx[2] = -1413.4;
  else{
    std::cout << "[ERROR]: " << "Incorrect injector position. Please enter B1-5." << std::endl;
    exit(0);
  }
  
  return vtx;
}


std::vector<double> GetTargetLocation(std::string inj){

  std::vector<double> vtx(3);

  if(inj == "B1"){      vtx[0] = -1595; vtx[1] = -550.6; vtx[2] = 1223;  }
  else if(inj == "B2"){ vtx[0] = -1567; vtx[1] = -620;   vtx[2] = 753;   }
  else if(inj == "B3"){ vtx[0] = -1575; vtx[1] = -553.6; vtx[2] = -215;  }
  else if(inj == "B4"){ vtx[0] = -1515; vtx[1] = -810;   vtx[2] = -756.2;}
  else if(inj == "B5"){ vtx[0] = -1500; vtx[1] = -600;   vtx[2] = -1413; }
  else{
    std::cout << "[ERROR]: " << "Incorrect injector position. Please enter B1-5." << std::endl;
    exit(0);
  }
  
  return vtx;
}


#endif
