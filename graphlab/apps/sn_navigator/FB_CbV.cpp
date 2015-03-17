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
    
    std::multimap<std::string, std::string>  att_values;

    single_message(): att_values() { }

   void save(graphlab::oarchive& oarc) const {
                size_t num = att_values.size();
                oarc << num;
                for(std::multimap<std::string, std::string>::const_iterator iter = att_values.begin();
                       iter != att_values.end(); ++iter){
                               oarc << iter->first;
                               oarc << iter->second;
                }

   }
   void load(graphlab::iarchive& iarc) {
          size_t num;
          std::string vk,v;
          iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk >> v;
                  att_values.insert(std::make_pair(vk,v));
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
        bool is_centroid;        
        std::vector<std::string> att_names;
        std::multimap<std::string, std::string>  att_values;
        std::multimap<std::string, std::vector<std::string> > att_v_n_freq;
        std::map<std::string, size_t> att_group;

     vdata() :
                att_names(), att_values(), att_v_n_freq(), att_group(), is_centroid(false) {}

      explicit vdata(std::vector< std::string> a_n, std::multimap<std::string, std::string> a_v,  bool isc) :
               att_names(a_n), att_values(a_v), att_v_n_freq(), att_group(), is_centroid(isc) {}

        void save(graphlab::oarchive& oarc) const {
               
                //store list of friends
                oarc << is_centroid;
                
                //CbIV
                size_t num = att_names.size();
                oarc << num;
                for(std::vector<std::string>::const_iterator iter = att_names.begin();
                       iter != att_names.end(); ++iter){
                               oarc << *iter;
                }
                
                num = att_values.size();
                oarc << num;
                for(std::multimap<std::string, std::string>::const_iterator iter = att_values.begin();
                       iter != att_values.end(); ++iter){
                               oarc << iter->first;
                               oarc << iter->second;
                       //std::cout <<"key: "<<iter->first<<", val: "<<iter->second<< "\n";
                }
                num = att_group.size();
                oarc << num;
                for(std::map<std::string, size_t>::const_iterator iter = att_group.begin();
                       iter != att_group.end(); ++iter){
                               oarc << iter->first;
                               oarc << iter->second;
                }
                num = att_v_n_freq.size();
                oarc << num;
                for(std::multimap<std::string, std::vector<std::string> >::const_iterator iter = att_v_n_freq.begin();
                       iter != att_v_n_freq.end(); ++iter){
                               oarc << iter->first;
                               std::vector<std::string> vlist =  iter->second;
                               for(size_t vi =0; vi<vlist.size(); vi++)
                                   oarc << vlist[vi];;
                }
                



        }
     void load(graphlab::iarchive& iarc) {
               
               std::string v,vk;
               size_t num = 0;

               iarc >> is_centroid;

               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk ;
                  att_names.push_back(vk);
               }

               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk >> v;
                  att_values.insert(std::make_pair(vk,v));
                  //if(vk != "interest")
                  //std::cout <<"key: "<<vk<<", val: "<<v<< "\n";
               }

               iarc >> num;
               size_t freq=0;
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk >> freq;
                  att_group.insert(std::make_pair(vk,freq));
               }

               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  std::vector<std::string> vs;
                  iarc >> vk ;
                  //Att1value
                  iarc >> v;
                  vs.push_back(v);
                  //Att2value
                  iarc >> v;
                  vs.push_back(v);
                  //Att3value
                  iarc >> v;
                  vs.push_back(v);

                  att_v_n_freq.insert(std::make_pair(vk,vs));
               }
        }

 };


///////////////////////////////global functions/////////////////////////////////
bool isCentroid(  graphlab::vertex_id_type vid ){
   bool chk = false;
   std::map< graphlab::vertex_id_type,  graphlab::vertex_id_type> centroids;
   std::map< graphlab::vertex_id_type,  graphlab::vertex_id_type>::iterator it;

   centroids[1019578116] = 1019578116;
   centroids[721547651] = 721547651;
   centroids[780219903] = 780219903;
   centroids[615038675] = 615038675;
   centroids[1034232372] = 1034232372;
   centroids[1448183152] = 1448183152;
   centroids[775333694] = 775333694;
   centroids[543733258] = 543733258;
   centroids[610148917] = 610148917;
   centroids[687892729] = 687892729;
   centroids[656597045] = 656597045;
   centroids[659229259] = 659229259;
   centroids[651308035] = 651308035;
   centroids[595587084] = 595587084;
   centroids[555825818] = 555825818;
   centroids[543902003] = 543902003;
   centroids[588672869] = 588672869;
   centroids[13611296] = 13611296;
   centroids[525476270] = 525476270;

   it=centroids.find(vid);
   if(it != centroids.end())
     chk = true;

   return chk;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//distributed graph definition
typedef graphlab::distributed_graph<vdata, edata> graph_type;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//#################################### Send data Items part #####################################
//

struct data_gatherer_directN {

        std::multimap<std::string, std::string> att_values;

        data_gatherer_directN () : att_values() {}
        void save(graphlab::oarchive& oarc) const {
            size_t num = att_values.size();
                oarc << num;
                for(std::multimap<std::string, std::string>::const_iterator iter = att_values.begin();
                       iter != att_values.end(); ++iter){
                               oarc << iter->first;
                               oarc << iter->second;
                      std::cout <<"key: "<<iter->first<<", val: "<<iter->second<< "\n";
                }

        }
        void load(graphlab::iarchive& iarc) {
             size_t num;
             std::string vk, v;
             att_values.clear();
             iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk >> v;
                  att_values.insert(std::make_pair(vk,v));
               }

        }
       data_gatherer_directN& operator+=(const data_gatherer_directN& other) {
              for(std::multimap<std::string, std::string>::const_iterator iter = other.att_values.begin();
                       iter != other.att_values.end(); ++iter){
                  att_values.insert(std::make_pair(iter->first,iter->second));
                   //if(iter->first!="interest")
                   //std::cout <<"key: "<<iter->first<<", val: "<<iter->second<< "\n";
              }
              return *this;
         }
};

class collect_attributes:public graphlab::ivertex_program<graph_type,data_gatherer_directN,data_message>, public graphlab::IS_POD_TYPE {

public:
 
       std::multimap<std::string, std::string> att_values;

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
          if(isCentroid (vertex.id())){
             /*for(std::multimap<std::string, std::string>::const_iterator iter = vertex.data().att_values.begin();
                       iter != vertex.data().att_values.end(); ++iter) {
                       att_n_v.att_values.insert(std::make_pair(iter->first,iter->second));
                       //if( iter->first =="education")
                       //std::cout <<"key: "<<iter->first<<", val: "<<iter->second<< "\n";
              }*/
              
             if(edge.target().id() == vertex.id() && (edge.source().data().att_values.size() >0))               
               for(std::multimap<std::string, std::string>::iterator iter = edge.source().data().att_values.begin(); 
                   iter!=edge.source().data().att_values.end(); iter++){
                   std::string vk = iter->first;
                   std::string v = iter->second;
                   att_n_v.att_values.insert(std::make_pair(vk, v));
                   //if( vk !="interest")
                   //std::cout <<"key: "<<vk<<", val: "<<v<< "\n";
              }
          }

        
          return att_n_v;
        }

        void scatter(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
        }


   //====================================================================================
   void apply(icontext_type& context, vertex_type& vertex,const gather_type& total) {
         
        //context.cout() << "apply it: " << context.iteration() << "\n" ;
        if(context.iteration() == 0 && isCentroid( vertex.id()) ){
           //first loop over gathered data and add them to vertex att_values
           //context.cout()<<"node id: "<<vertex.id()<<", gather size: "<< total.att_values.size() << ", node att size: " << vertex.data().att_values.size() <<   "\n" ;
           
           for(std::multimap<std::string, std::string>::const_iterator iter = total.att_values.begin(); iter!= total.att_values.end(); ++iter){
               vertex.data().att_values.insert(std::make_pair(iter->first, iter->second));
               //context.cout()<<"key: "<<iter->first<<", val: "<<iter->second<< "\n";
           }
          
           //calculting the frequencies using node's att_values
           for(size_t i=0;i<vertex.data().att_names.size(); i++){
               std::string a_n1 = vertex.data().att_names[i];
               std::vector<std::string> triple;
               for(size_t j=1; j<vertex.data().att_names.size(); j++){
                  std::string a_n2 = vertex.data().att_names[j];   
                  if(a_n1 != a_n2){
                  context.cout() <<"an1: "<<a_n1 << " with an2: " << a_n2 << "\n";
                  std::pair<std::multimap<std::string, std::string>::iterator,std::multimap<std::string, std::string>::iterator > 
                       iter_n1 = vertex.data().att_values.equal_range(a_n1);

                  std::pair<std::multimap<std::string, std::string>::iterator,std::multimap<std::string, std::string>::iterator >
                       iter_n2 = vertex.data().att_values.equal_range(a_n2);

                  std::vector<std::string> val_n1, val_n2;
                  for (std::multimap<std::string, std::string>::iterator it=iter_n1.first; it!=iter_n1.second; ++it){
                      val_n1.push_back(it->second);
                      context.cout() <<"values for an1: " << it->second << "\n";
                  }
                  for (std::multimap<std::string, std::string>::iterator it=iter_n2.first; it!=iter_n2.second; ++it){
                      val_n2.push_back(it->second);
                      context.cout() <<"values for an2: " << it->second << "\n";
                  }
                  std::string value1, value2;
                  size_t freq = 1;
                  size_t len = val_n1.size();
                  if(len > val_n2.size())
                     len = val_n2.size();
                  for(size_t vali =0; vali<len; vali ++){
                     value1 = val_n1[vali];
                     value2 = val_n2[vali];
                     //context.cout() <<"for freq val1: "<<value1 << " with val2: " << value2 << "\n";        
                     //start to find repetition and increase counter
                     //check if counting is done before for same values
                     std::multimap<std::string,std::vector<std::string> >::iterator fit = vertex.data().att_v_n_freq.find(value1+ "&" + value2);
                     if(fit == vertex.data().att_v_n_freq.end()){
                        //context.cout() <<"new combination for freq \n";
                        triple.clear();
                        triple.push_back(a_n1);
                        triple.push_back(a_n2);     
                        for(size_t index = vali +1; index<len; index++)
                           if(value1 == val_n1[index] && value2 == val_n2[index])
                              freq ++; 
                        std::stringstream convert;
                        convert << freq;
                        triple.push_back(convert.str());
                        vertex.data().att_v_n_freq.insert(std::make_pair(value1+ "&" + value2, triple));
                     }  
                  }
                  }
               }
           }
        }
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
    std::vector<graphlab::vertex_id_type> frnds;
    std::string attribute, token, token_sub, token_list;
    bool isCent = false;
    std::vector<std::string> att_names;
    std::multimap<std::string, std::string>  att_values;

    att_names.push_back("gndr");
    att_names.push_back("fname");
    att_names.push_back("education");
    //att_names.push_back("edu_org");
    //att_names.push_back("edu_city");
    //att_names.push_back("edu_cntry");
    att_names.push_back("job");
    att_names.push_back("employer");
    //att_names.push_back("job_city");
    //att_names.push_back("job_cntry");
    att_names.push_back("interest");
    //att_names.push_back("home_city");
    att_names.push_back("home_cntry");
    //att_names.push_back("currloc_city");
    att_names.push_back("currloc_cntry");

/////////////////////////////////////////////////Load data    
    strm >> vid;
    strm >> attribute;
   
    //graph.add_vertex(vid);
    isCent = isCentroid (vid);
    if(isCent) {
       std::istringstream ss(attribute);
       while(std::getline(ss, token, ',')) {
          size_t fID = atoi(token.c_str());
          if(fID != vid){
             frnds.push_back(fID);
             //graph.add_edge(fID, vid);
          }
       }
    }
 
   strm >> attribute;
   std::istringstream ss1(attribute);
   while(std::getline(ss1, token, ','))
      att_values.insert(std::make_pair("job",token));
 
   strm >> attribute;
   std::istringstream ss2(attribute);
   while(std::getline(ss2, token, ','))
       att_values.insert(std::make_pair("employer",token));
   
   //strm >> attribute;
   //std::istringstream ss9(attribute);
   //while(std::getline(ss9, token, ','))
      //att_values.insert(std::make_pair("job_city",token));

   //strm >> attribute;
   //std::istringstream ss3(attribute);
   //while(std::getline(ss3, token, ','))
      //att_values.insert(std::make_pair("job_cntry",token));

   strm >> attribute;
   std::istringstream ss4(attribute);
   while(std::getline(ss4, token, ',')){
      int index =0;
      std::string cat;
      std::istringstream ss12(token);
       while(std::getline(ss12, token_sub, '$')){
         if(index ==0){
            cat = token_sub;
            index ++;
         }
         else{
          att_values.insert(std::make_pair("education",cat + "_" + token_sub));
          //std::cout << "edu " << cat + "_" + token_sub << "\n"; 
         }
       }
   }


   //strm >> attribute;
   //std::istringstream ss5(attribute);
   //while(std::getline(ss5, token, ','))
      //att_values.insert(std::make_pair("edu_org",token));

   //strm >> attribute;
   //std::istringstream ss6(attribute);
   //while(std::getline(ss6, token, ','))
       //att_values.insert(std::make_pair("edu_city",token));

   //strm >> attribute;
   //std::istringstream ss7(attribute);
   //while(std::getline(ss7, token, ','))
      //att_values.insert(std::make_pair("edu_cntry",token));
   
   strm >> attribute;
   std::istringstream ss8(attribute);
   //#cat1@tok1,tok2#cat2@tok21,tok22
   //split categories
   while(std::getline(ss8, token, '#')){
      int index =0;
      std::string cat;
      std::istringstream ss10(token);
      while(std::getline(ss10, token_sub, '@')){
         if(index ==0){
            cat = token_sub;
            index ++;
         }
         else{
             std::istringstream ss11(token_sub);
             while(std::getline(ss11, token_list, ',')){   
                 att_values.insert(std::make_pair("interest",cat + "_" + token_list));  
                 //std::cout << "interest " << cat + "_" + token_list << "\n";
             }
         }
      }
   }
   
   strm >> attribute;
   //std::cout << "fname " << attribute << "\n";
   att_values.insert(std::make_pair("fname",attribute));

   strm >> attribute;
   //std::cout << "gndr " << attribute << "\n";
   att_values.insert(std::make_pair("gndr",attribute));

   //strm >> attribute;
   //att_values.insert(std::make_pair("home_city",attribute));
 
   strm >> attribute;
   //std::cout << "hc " << attribute << "\n";
   att_values.insert(std::make_pair("home_cntry",attribute));

   //strm >> attribute;
   //att_values.insert(std::make_pair("currloc_city",attribute));

   strm >> attribute;
   //std::cout << "curc " << attribute << "\n";
   att_values.insert(std::make_pair("currloc_cntry",attribute));
   
  //for(std::multimap<std::string, std::string>::iterator iter = att_values.begin(); iter!=att_values.end(); iter++)
           //std::cout << "key: " << iter->first << ", val: " << iter->second << "\n"; 


///////////////////////////////////////////Add node to the graph
    graph.add_vertex(vid, vdata(att_names,att_values,isCent));
    if(isCent)
      for(int  i=0; i< frnds.size(); i++)
         graph.add_edge(frnds[i], vid);
   


   return true;
 }


//#################################  Graph Writer ######################################//
//#################################  Graph Writer ######################################//
struct graph_writer {
  std::string save_vertex(graph_type::vertex_type v) {


        std::stringstream strm;
        if(isCentroid(v.id())){
        strm <<  v.id() << "\n";
        for(std::multimap<std::string, std::vector<std::string> >::iterator iter = v.data().att_v_n_freq.begin();
                       iter != v.data().att_v_n_freq.end(); ++iter){
                               strm << iter->first << "\t";
                               std::vector<std::string> vlist =  iter->second;
                               for(size_t vi =0; vi<vlist.size(); vi++)
                                    strm << vlist[vi] << "\t";
                                strm <<"\n";
                }
       } 

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

   
   size_t MAX_ITERATIONS=1;
   clopts.get_engine_args().set_option("max_iterations", MAX_ITERATIONS);
   graphlab::omni_engine<collect_attributes> engine2(dc, graph, "sync", clopts);
   engine2.signal_all();
   engine2.start();

   graph.save("FB_Data.txt", graph_writer(),
             false,    // do not gzip
             true,     // save vertices
             false);   // do not save edges
   
   graphlab::mpi_tools::finalize();



   return EXIT_SUCCESS;
  }



