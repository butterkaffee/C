//Author: Michael Fitzke

#include <vector>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath> 

using namespace std;

struct Connection
{
  double weight;
  double deltaWeight; 
}  

typedef vector<Neuron> Layer; 
//*****************CLASS NEURON****************

class Neuron
{
public:
  Neuron(unsigned numOutputs, unsigned myIndex);
  void setOutputVal(double val){m_outputVal = val;}
  double getOutputVal(void) const {return n_outputVal;}
  void feedForward(const Layer &prevLayer);
 
private;
  static double transferFunction(double x);
  static double transferFunctionDerivative(double x); 
  static double randomWeight(void) {return rand() /double(RAND_MAX)}; 
  double m_outputVal;
  vector<Connection> m_outputWeights;
  unsigned m_myIndex;
}

Neuron::transferFunction(double x)
{
  //tanh Funktion

  return tanh(x);
}

double Neuron::transferFunctionDerivative(double x)
{
  return 1 - tanh(x)*tanh(x);
}

Neuron::feedForward(const Layer &prevLayer){
  double sum = 0.0;


  for (unsigned n=0; n<prevLayer.size(); ++n){
     sum += prevLayer[n].getOutputVal()* 
             prevLayer[n].m_outputWeights[m_myIndex].weight;
   }
   
   m_outputVal = Neuron::transferFunction(sum); 

}


Neuron::Neuron(unsigned numOutputs, unsigned myIndex){
  for (unsigned c = 0; c < numOutputs, ++c){
     m_outputWeights.push_back(Connection() );
     m_outputWeights.back().weight = randomWeight();     
  }
  m_myIndex = myIndex; 
}

//******************CLASS NET***************
class Net
{
public:
  Net(const vector<unsigned> &topology);
  void feedForward(const vector<double> &inputVals);
  void backProp(const vector<double> &targetVals); 
  void getResults(vector<double> &resultVals) const {}; 


private: 
  vector<Layer> m_layers; //m_layers[layerNum][neuronNum]
  double m_error;
  double m_recentAverageError; 
  double m_recentAverageSmoothingFactor;  
};

Net::void backProp(const vector<double> &targetVals)
{
  //Calc overall net Error

  Layer &outputLayer = m_layers.back(); 
  m_error = 0.0;
  
  for (unsigned n= 0; n<outpuLayer.size() -1; ++n){
    double delta = targetVals[n] - outputLayer[n].getOutputVal();
    m_error += delta*delta;
  }
  m_error /= outputLayer.size() -1;
  m_error = sqrt(m_error); 
  
  //Implement a recent average measurement
  
  m_recentAverageError = 
    (m_recentAverageError * m_recentAverageSmoothingFactor + m_error)
    / (m_recentAverageSmoothingFactor +1.0); 

  //Calc output layer Gradients 
  
  for (unsigned n = 0; n < outputLayer.size()-1, ++n){
     outputLayer[n].calcOutputGradients(targetVals[n]); 

  }   

  //Calc gradients on hidden Layer
  
  for (unsigned layerNum = m_layers.size() - 2; layerNum > 0; --layerNum){
  
    Layer &hiddenLayer = m_layers[layerNum]; 
    Layer &nextLayer = m_layer[layerNum +1];
  
    for(unsigned n = 0; n < hiddenLayer.size(); ++n) {
    
      hiddenLayer[n].calcHiddenGradients(nextLayer);

    }
  
  for (unsigned layerNum = m_layers.size() - 1; layerNum >0; --layerNum){

  }
  }

   


  //For all layers from outputs to first hidden layer 
  // update connection weights
} 
Net::feedForward(const vector<double> &inputVals)
{
   assert(inputVals.size() == m.layers[0].size() - 1 );
   

   //Input Neurons setzen 
   for (unsigned i = 0; i < inputVals.size(); ++i){
     m_layers[0][i].setOutputVal(inputVals[i]);
  }

      //Forward Porpagate
     for (unsigned layerNum = 1; layerNum < m_layers.size(); ++layerNum){
        Layer &prevLayer = m_layers[layerNum - 1];
        for (unsigned n = 0; n<m_layers[layerNum].size() - 1; ++n) {
           //neuron Ebene
           m_layers[layerNum][n].feedForward(prevLayer);
    }
  }
}



Net::Net(const vector<unsigned> &topology)
{
  unsigned numLayers = topology.size();

  for (unsigned layerNum = 0; layerNum < numLayers; ++layerNum){
     m_layers.push_back(Layer()); 
     unsigned numOutputs = layerNum == topology.size() - 1 ? 0: topology[layerNum + 1];

    // we have new layer, now fill it with ith neurons and add bias

     for (unsigned neuronNum = 0; neuronNum <= topology[layerNum]; ++neuronNum){
    
     // CONTAINER<X>.back() gives back most recently added object 
        m_layers.back().push_back(Neuron(numOutputs, neuronNum));
        cout << "Made a Neuron" << endl;  
    }
  }
}

int main()
{
  // e.g. {3,2,1}

  vector<unsigned> topology;
  topology.push_back(3);
  topology.push_back(2);
  topology.push_back(1); 
  Net myNet(topology);

  vector<double> inputVals;
  myNet.feedForward(inputVals);

  vector<double> targetVals; 
  myNet.backProp(targetVals);

  vector<double> resultVals; 
  myNet.getResults(resultVals);

}
