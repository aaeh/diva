/* This software is distributed under the MIT license (see LICENSE file in the
 * root directory)
 *
 * Copyright (c) 2013 Amira A. Soliman
 */

#include <graphlab.hpp>
#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/ui/metrics_server.hpp>
#include <graphlab/macros_def.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <algorithm>    // std::reverse
/*
*#################################### Data Structure Part #####################################
*/

size_t VIEW_SIZE = 50;
size_t PAYLOAD_SIZE = 10;


//////////////////////////////////////////////////////////
//gossip message exchanged between social graph nodes
//////////////////////////////////////////////////////////

struct single_message {
    graphlab::vertex_id_type sender_id;
    graphlab::vertex_id_type reciever_id;  
    std::vector<graphlab::vertex_id_type> path;

    single_message(): path() { }

   void save(graphlab::oarchive& oarc) const {
                oarc <<sender_id<<reciever_id;
                size_t num = path.size();
                oarc << num;
                for(std::vector<graphlab::vertex_id_type>::const_iterator iter = path.begin();
                       iter != path.end(); ++iter){
                               oarc << *iter;
                }




   }
   void load(graphlab::iarchive& iarc) {
          size_t num=0;
          iarc>>sender_id>>reciever_id;
          graphlab::vertex_id_type f;
          iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> f;
                  path.push_back(f);
               }
   }
};

// message object exchanged using ivertex program
struct data_message: public graphlab::IS_POD_TYPE {
  std::vector<single_message> msgs;
  double value;

  data_message():  value(1.0), msgs() { }

  //
  // returns a priority to the reception of the message (used by the asynchronous engines)
  //
  double priority() const { return std::fabs(value); }

  //
  // allow message merging across the network
  // GraphLab always merges messages destined to the same vertex.
  //
  data_message& operator+=(const data_message& other) {

          for(size_t i = 0;i < other.msgs.size(); ++i){
               msgs.push_back(other.msgs[i]);
          }

          return *this;
  }
  void save(graphlab::oarchive& oarc) const {

        size_t num_msgs = msgs.size();
        oarc << num_msgs;
        for(size_t k = 0;k < num_msgs; ++k){
                 msgs[k].save(oarc);
        }
  }

  void load(graphlab::iarchive& iarc) {

          size_t num_msgs = 0;
          msgs.clear();

          iarc >> num_msgs;

          for(size_t k = 0;k < num_msgs; ++k){
                  single_message msg;
                  msg.load(iarc);
                  msgs.push_back(msg);
          }
  }
};

//////////////////////////////////////////////////////////
//edge
struct edata {
  graphlab::edge_id_type weight;

  edata():weight(0) {
  }

  explicit edata(graphlab::edge_id_type w): weight(w) {
  }

  void save(graphlab::oarchive& oarc) const {
    oarc << weight;
  }

  void load(graphlab::iarchive& iarc) {
    iarc >> weight;
  }

};
////////I//////////////////////////////////////////////////
//vertex
struct vdata {

        // list of direct friends
                
        size_t commID;
        std::vector<graphlab::vertex_id_type> friends;
        

     vdata() : friends(), commID(0) {}
      explicit vdata(std::vector<graphlab::vertex_id_type> f,  size_t cid) :
          friends(f), commID(cid) {}


        void save(graphlab::oarchive& oarc) const {
              
                //store list of friends
                
                oarc << commID;
                size_t num = friends.size();
                oarc << num;
                for(size_t i = 0;i < num; ++i){
                        oarc << friends[i];
                }


             
        }
     void load(graphlab::iarchive& iarc) {
               
               size_t num = 0, f=0;
               
             
               iarc >> commID;

               iarc >> num;
                for(size_t i = 0; i < num; i++) {
                  graphlab::vertex_id_type f_id = 0;
                  iarc >> f_id;
                  friends.push_back(f_id);
                }


        }
 };


///////////////////////////////global functions/////////////////////////////////
graphlab::vertex_id_type find_max(std::map<graphlab::vertex_id_type, size_t> vec, size_t maxC){
    graphlab::vertex_id_type max=maxC;
    
    for(std::map<graphlab::vertex_id_type, size_t>::iterator it=vec.begin(); it != vec.end(); it++){
       if(it->second >  max){
          max =it->second;
       }
    }
    return max;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//distributed graph definition
typedef graphlab::distributed_graph<vdata, edata> graph_type;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//#################################### Send data Items part #####################################
//

struct data_gatherer_directN {

        std::map<graphlab::vertex_id_type, size_t> att_values;


       data_gatherer_directN () : att_values() {}
        void save(graphlab::oarchive& oarc) const {
           size_t num = att_values.size();
                oarc << num;
                for(std::map<graphlab::vertex_id_type, size_t>::const_iterator iter = att_values.begin();
                       iter != att_values.end(); ++iter){
                               oarc << iter->first;
                               oarc << iter->second;
                      
                }
            
        }
        void load(graphlab::iarchive& iarc) {
             size_t num, c;
              att_values.clear();
             graphlab::vertex_id_type f;
            
             iarc >> num;
             for(int j=0; j<num;j++){
               iarc >> f>> c;
               att_values.insert(std::make_pair(f,c));
             }
            }
   
       data_gatherer_directN& operator+=(const data_gatherer_directN& other) {
              for(std::map<graphlab::vertex_id_type, size_t> ::const_iterator iter = other.att_values.begin();
                       iter != other.att_values.end(); ++iter){
                    att_values.insert(std::make_pair( iter->first,iter->second));
                   //if(iter->first!="interest")
                   //std::cout <<"key: "<<iter->first<<", val: "<<iter->second<< "\n";
              }

              return *this;
         }

};

class collect_attributes:public graphlab::ivertex_program<graph_type,data_gatherer_directN,data_message>, public graphlab::IS_POD_TYPE {

public:
 
       std::map<graphlab::vertex_id_type, size_t> att_values;
       

        edge_dir_type gather_edges(icontext_type& context,
                              const vertex_type& vertex) const {
                return graphlab::ALL_EDGES;
        }

        edge_dir_type  scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const {
                return graphlab::NO_EDGES;
        }

        data_gatherer_directN gather(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
          data_gatherer_directN att_n_v;
             /*for(std::multimap<std::string, std::string>::const_iterator iter = vertex.data().att_values.begin();
                       iter != vertex.data().att_values.end(); ++iter) {
                       att_n_v.att_values.insert(std::make_pair(iter->first,iter->second));
                       //if( iter->first =="education")
                       //std::cout <<"key: "<<iter->first<<", val: "<<iter->second<< "\n";
              }*/
              att_n_v.att_values.insert(std::make_pair( edge.source().id() , edge.source().data().commID));
          return att_n_v;
        }

        void scatter(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
        }


   //====================================================================================

   void init(icontext_type& context,
                const vertex_type& vertex, const data_message &msg) {
      if(vertex.data().friends.size() >0)    
         if(context.iteration() == 0 ){
             
          }
      } 
   }


   void apply(icontext_type& context, vertex_type& vertex,const gather_type& total) {
            //aggregate community results
 }
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//############################ line parser to load data from file ##################################################//

bool line_parser(graph_type& graph,
                   const std::string& filename,
                   const std::string& textline)
 {
    if(textline.empty()) return true;

    std::stringstream strm(textline);

    graphlab::vertex_id_type vid;
    std::vector<graphlab::vertex_id_type> frnds, c_frnds;
    std::string attribute, token, token_sub, token_list;
    bool  hasFrnds = false, hasAttrib=false;
    size_t cid;
    std::vector<std::string> att_names;
    std::multimap<std::string, std::string>  att_values,c_att_values;

    att_names.push_back("school");
    att_names.push_back("major");
    att_names.push_back("employer");
    att_names.push_back("places_lived");

/////////////////////////////////////////////////Load data    
    strm >> vid;
    strm >> cid;
    
    strm >> attribute;
    if(attribute != "null"){
       hasFrnds = true;
       std::istringstream ss(attribute);
       while(std::getline(ss, token, '\\')) {
          size_t fID = atoi(token.c_str());
          if(fID != vid){
             frnds.push_back(fID);
          }
       }
    }
   

///////////////////////////////////////////Add node to the graph
  //if(cid == 61935 || vid == 61935 || vid==4851 || vid==438 || vid==698){
   graph.add_vertex(vid, vdata(frnds,cid));
    //if(isCent)
      for(int  i=0; i< frnds.size(); i++)
         graph.add_edge(frnds[i], vid);


   return true;
 }


//#################################  Graph Writer ######################################//
//#################################  Graph Writer ######################################//
struct graph_writer {
  std::string save_vertex(graph_type::vertex_type v) {


        std::stringstream strm;
        std::vector<std::string> att_vec;
        std::vector<size_t> att_freq;
        if(v.data().friends.size()>0){    
          strm << "============================================== \n";
          strm <<  v.id() << " , " << v.data().friends.size() << "\n";
          
          /*for(size_t i=0; i< v.data().att_group_clstr.size(); i++){
             double supp = double(v.data().freq_clstr[i])/double(v.data().supp_clstr[i]);
             double conf = supp/double(v.data().conf_clstr[i]);
             strm << v.data().att_group_clstr[i] << " , " << v.data().freq_clstr[i] << " , " << supp << " , " << conf << "\n";
          }*/
       //}
       //for(size_t i=0; i< v.data().att_group.size(); i++)
            // strm << v.data().att_group[i] << " , " << v.data().freq[i] << "\n";

    return strm.str();
  }
  std::string save_edge(graph_type::edge_type e) { return ""; }
};

//######################################## main function #########################################//
//######################################## main function #########################################//
  int main(int argc, char** argv) {

        graphlab::command_line_options clopts("Dosn: fill neighbors map");
        std::string file_name;
        clopts.attach_option("graph", file_name, "Graph input. reads all graphs matching prefix*");

        if(!clopts.parse(argc, argv)) return EXIT_FAILURE;
        if (file_name == "") {
                std::cout << "--graph is not optional\n";
                clopts.print_description();
                return EXIT_FAILURE;
        }

        // Initialize control plane using mpi
        graphlab::mpi_tools::init(argc, argv);
        graphlab::distributed_control dc;
        graphlab::launch_metric_server();

        // load graph
        graph_type graph(dc, clopts);      //graph_type graph(dc);
         
        graph.load(file_name, line_parser); //load_from_hdfs("hdfs://graphLabLarge1:154310/data/graph.txt", line_parser);
        //graph.load_format(file_name, "tsv");
        graph.finalize();
        dc.cout() << "Number of vertices: " << graph.num_vertices() << std::endl
              << "Number of edges:    " << graph.num_edges() << std::endl;

   
   size_t MAX_ITERATIONS=9;
   clopts.get_engine_args().set_option("max_iterations", MAX_ITERATIONS);
   graphlab::omni_engine<collect_attributes> engine2(dc, graph, "sync", clopts);
   engine2.signal_all();
   engine2.start();
   std::cout << "Runtime: " << engine2.elapsed_seconds();
   graph.save(file_name, graph_writer(),
             false,    // do not gzip
             true,     // save vertices
             false);   // do not save edges
  
   graphlab::mpi_tools::finalize();



   return EXIT_SUCCESS;
  }



