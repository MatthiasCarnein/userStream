#include <Rcpp.h>
#include <limits>

#define VERBOSE 0

using namespace Rcpp;


class CF {
public:
  Rcpp::NumericVector LS;
  Rcpp::NumericVector SS;
  double weight;
  int t;

  // constructor
  CF(Rcpp::NumericVector LS, Rcpp::NumericVector SS, double weight, int t){
    this->LS = LS;
    this->SS = SS;
    this->weight = weight;
    this->t = t;
  }

  // copy constructor
  CF(CF* cf){
    this->LS = Rcpp::clone(cf->LS);
    this->SS = Rcpp::clone(cf->SS);
    this->weight = cf->weight;
    this->t = cf->t;
  }


  double getWeight(){
    return this->weight;
  }

  double distance(CF* cf){
    return sqrt(sum(pow(this->getCentroid() - cf->getCentroid(), 2)));
  }

  Rcpp::NumericVector getCentroid(){
    return Rcpp::clone(this->LS) / this->weight;
  }

  double getRadius() {
    return sum(this->SS/this->weight - pow((this->LS/this->weight),2)); // variance
  }

  void merge(CF* cf, int t, double lambda) {
    // fade (with updating time)
    this->fade(t, lambda);
    cf->fade(t, lambda);

    // sum components
    this->weight = this->weight+cf->weight;
    this->LS = this->LS + cf->LS;
    this->SS = this->SS + cf->SS;
  }


  void remove(CF* cf, int t, double lambda){
    this->fade(t, lambda);
    cf->fade(t, lambda);

    // subtract components
    this->weight = this->weight - cf->weight;
    this->LS = this->LS - cf->LS;
    this->SS = this->SS - cf->SS;

  }

  void scale(Rcpp::NumericVector mean, Rcpp::NumericVector SD){
    this->LS = this->LS / SD;
    this->SS = this->SS / pow(SD,2);
  }

  void unscale(Rcpp::NumericVector mean, Rcpp::NumericVector SD){
    this->LS = this->LS * SD;
    this->SS = this->SS * pow(SD,2);
  }

  void fade(int t, double lambda) {
    if(lambda == 0) return;

    // apply fading
    this->weight = this->weight * pow(2, (-lambda * (t-this->t)));
    this->LS = this->LS * pow(2, (-lambda * (t-this->t)));
    this->SS = this->SS * pow(2, (-lambda * (t-this->t)));

    this->t = t;
  }

};





class UserStream {
public:

  double r;
  double lambda;
  int tnow;
  int tgap;
  double omega;
  int d;
  int n;
  std::vector<CF*> micro;
  std::map<int, int> assignment;
  std::map<int, CF*> lastCFs;

  // keep track of mean and SD for scaling
  Rcpp::NumericVector mean;
  Rcpp::NumericVector dev;
  int count;

  Rcpp::NumericVector currentMean;
  Rcpp::NumericVector currentSD;

  int numScales;


  UserStream(double r, double lambda, int tgap){
    this->r = r;
    this->lambda = lambda;
    this->tnow = 0;
    this->count=0;
    this->tgap = tgap;
    this->omega = 0;
    if(lambda !=0) this->omega = pow(2, (-1*lambda * tgap));
    this->d = 0;
    this->n = 0;
    this->numScales =1;
  }




  Rcpp::IntegerMatrix getAssignment(){
    this->fadeAll();

    Rcpp::IntegerMatrix mapResult(this->assignment.size(),2);
    int i=0;
    for (std::map<int, int>::iterator it = this->assignment.begin(); it != this->assignment.end(); it++ ){
      mapResult(i,0)=it->first;
      mapResult(i,1)=it->second+1;
      i++;
    }
    return mapResult;
  }

  Rcpp::NumericMatrix getLastCFs(){
    Rcpp::NumericMatrix mapResult(this->lastCFs.size(), this->d);
    int i=0;
    for (std::map<int, CF*>::iterator it = this->lastCFs.begin(); it != this->lastCFs.end(); it++ ){
      mapResult(i,_)=it->second->getCentroid();
      i++;
    }
    return mapResult;
  }


  Rcpp::NumericMatrix get_microclusters(bool unscale){
    this->fadeAll();

    Rcpp::NumericMatrix x(this->micro.size(),d);
    for(unsigned int i=0; i<this->micro.size();i++){

      CF* cf = new CF(this->micro[i]);
      if(unscale) cf->unscale(this->currentMean, this->currentSD);
      Rcpp::NumericVector centroid = cf->getCentroid();
      for(int j=0; j<centroid.size();j++){
        x(i,j)=centroid[j];
      }
      delete cf;
    }
    return(x);
  }


  Rcpp::NumericVector get_microweights(){
    this->fadeAll();

    Rcpp::NumericVector x(this->micro.size());
    for(unsigned int i=0; i<this->micro.size();i++){
      x[i]=this->micro[i]->getWeight();
    }
    return(x);
  }

  Rcpp::NumericVector getMean(){
    return Rcpp::clone(this->mean);
  }


  void fadeAll(){
    // fade Clusters
    for(int i=this->micro.size()-1; i>=0; i--){
      // Rcpp::Rcout << "Fade cluster " << i << std::endl;

      this->micro[i]->fade(this->tnow, this->lambda);
      this->removeInsufficient(i);
    }
  }




  Rcpp::NumericVector getVariance(){
    Rcpp::NumericVector var = this->dev/(this->count-1);
    for(int i = 0; i < var.size(); i++){
      if(var[i]<=0) var[i] = 1; // fix when no variance
    }
    return var;
  }


  Rcpp::NumericVector getSD(){
    return sqrt(this->getVariance());
  }


  void update(Rcpp::NumericVector newUser, int userId, int time) {

    // initialize here since dimension was unknown beforehand
    if(this->d==0){
      this->initialize(newUser.size());
    }

    this->n++;
    if(this->n % 1000 == 0) Rcpp::Rcout << n << std::endl;

    this->tnow = time;

    // check if user already known
    std::map<int, int>::iterator it = this->assignment.find(userId);
    if(it != this->assignment.end() && it->second >= 0){
      this->removeUser(userId, time); // remove user from its cluster
    }

    // insert updated user again
    this->insertUser(newUser, userId, time);
  }

  void removeUser(int userId, double time){
#if VERBOSE >= 1
    Rcpp::Rcout << userId << " - Remove User from Cluster " << this->assignment[userId]+1 << std::endl;
#endif

    CF* cf = this->lastCFs[userId]; // get last CF of user

    this->subtractFromScale(cf->getCentroid()); // adapt scaling parameters
    cf->scale(this->currentMean, this->currentSD); // scale it
    cf->fade(time, this->lambda);

    this->micro[this->assignment[userId]]->remove(cf, time, this->lambda); // remove from CF

    this->removeInsufficient(this->assignment[userId]); // check if CF can be removed
  }



  Rcpp::IntegerVector getClosest(Rcpp::NumericMatrix data){
    Rcpp::IntegerVector assignment(data.nrow());

    for(int i=0; i < data.nrow(); i++){
      Rcpp::NumericVector newUser = data(i,_);

      CF* cf = new CF(newUser, pow(newUser,2), 1, this->tnow); // create new CF for user
      cf->scale(this->currentMean, this->currentSD);

      int j = -1;
      double dist = std::numeric_limits<double>::infinity();
      for(unsigned int k=0; k < this->micro.size(); k++){
        double d = cf->distance(this->micro[k]);
        if(d < dist){
          dist = d;
          j=k;
        }
      }
      assignment(i) = j+1; // from C to R indexing
      delete cf;
    }
    return assignment;
  }



  void insertUser(Rcpp::NumericVector newUser, int userId, int time){
    this->addToScale(newUser); // adapt scaling parameters
    CF* cf = new CF(newUser, pow(newUser,2), 1, time); // create new CF for user

    this->lastCFs[userId] = new CF(cf); // remember CF (copy)

    // regularly adapt the scaling to the change in variance
    if(this->n % (int)pow(2,this->numScales)==0){
      if(this->numScales< 10) this->numScales++;
      this->adjustScale();
    }

    cf->scale(this->currentMean, this->currentSD); // scale it


    // search closest CF
    int assignedCluster;
    int j;
    int validCandidate=0;
    double dist;
    while(validCandidate==0){

      // if no micro cluster exists: insert
      if(!this->micro.size()){
#if VERBOSE >= 1
        Rcpp::Rcout << userId << " - Use User to create new Cluster " << this->micro.size()+1 << std::endl;
#endif
        this->micro.push_back(cf);
        // remember new user cluster
        this->assignment[userId] = this->micro.size()-1;;
        return;
      }

      // find closest CF
      j = -1;
      dist = std::numeric_limits<double>::infinity();
      for(unsigned int i=0; i < this->micro.size(); i++){
        double d = cf->distance(this->micro[i]);
        if(d < dist){
          dist = d;
          j=i;
        }
      }
      if(j==-1) Rcpp::Rcout << "Search for Cluster failed" << std::endl;
      // fade candidate and repeat if fading removed cluster
      this->micro[j]->fade(time, this->lambda);
      validCandidate = !this->removeInsufficient(j);
    }

    // copy mc to leave original intact
    CF* temp = new CF(this->micro[j]);

    // temporarily merge
    temp->merge(cf, time, this->lambda);

    // Check if radius has grown beyond threshold
    if(temp->getRadius() <= r){
#if VERBOSE >= 1
      Rcpp::Rcout << userId << " - Merge User into Cluster " << j+1 << " at dist " << dist << ". New Radius: " << temp.getRadius()<< " and weight: " << temp.getWeight()  << std::endl;
#endif
      this->micro[j] = temp;
      assignedCluster = j;
      delete cf;
    } else{
#if VERBOSE >= 1
      Rcpp::Rcout << userId << " - Use User to create new Cluster " << this->micro.size()+1 << std::endl;
#endif
      this->micro.push_back(cf);
      assignedCluster = this->micro.size()-1;
      delete temp;
    }
    // remember new user cluster
    this->assignment[userId] = assignedCluster;
  }



  // Welford’s method for computing variance
  void addToScale(Rcpp::NumericVector x){
    this->count++;

    // for every dimension
    for(int d=0; d<x.size(); d++){
      double mu2 = this->mean(d) + ((x(d) - this->mean(d))/this->count);
      this->dev(d) = this->dev(d) + ((x(d) - this->mean(d)) * (x(d) - mu2));
      this->mean(d) = mu2;
    }
  }

  void subtractFromScale(Rcpp::NumericVector x){
    this->count--;

    // for every dimension
    for(int d=0; d<x.size(); d++){
      double mu2 = mean(d) - ((x(d) - this->mean(d))/this->count);
      this->dev(d) = this->dev(d) - ((x(d) - this->mean(d)) * (x(d) - mu2));
      this->mean(d) = mu2;
    }

    // in case we removed everything, we start from scratch
    if(this->count==0){
      this->initialize(x.size());
    }
  }

  void adjustScale(){
    if(this->count<=1) return;

    Rcpp::NumericVector newMean = this->getMean();
    Rcpp::NumericVector newSD = this->getSD();

#if VERBOSE >= 1
    Rcpp::Rcout << "Adjust Scaling" << std::endl;
#endif

    for(int i=this->micro.size()-1; i>=0; i--){
      this->micro[i]->unscale(this->currentMean, this->currentSD);
      this->micro[i]->scale(newMean, newSD);
    }


    this->currentMean = newMean;
    this->currentSD = newSD;
  }

  void initialize(int d){
#if VERBOSE >= 1
    Rcpp::Rcout << "Initialize" << std::endl;
#endif
    this->d = d;
    this->mean = Rcpp::NumericVector(d);
    this->dev = Rcpp::NumericVector(d);
    this->currentMean = Rcpp::NumericVector(d);
    this->currentSD = Rcpp::NumericVector(d);
    for(int i=0; i<this->currentSD.size(); i++) this->currentSD[i]=1;
  }




  int removeInsufficient(int i){
// remove insufficient weight
    double weight = this->micro[i]->getWeight();
    if(weight <= this->omega){


      this->micro.erase(this->micro.begin()+i);

      // cluster assignment map
      std::map<int, int>::iterator it = this->assignment.begin();
      while(it != this->assignment.end()){
        if(it->second == i){
          int userId = it->first;
          std::map<int, CF*>::iterator userIt = this->lastCFs.find(userId);
          if(userIt != this->lastCFs.end()){
            // delete userIt->second;
            this->lastCFs.erase(userIt);
          }
          it->second = -1; // indicator for not assigned
        } else if(it->second > i){
          it->second--;
          it++;
        } else{
          it++;
        }
      }

#if VERBOSE >= 1
      Rcpp::Rcout << "Removed cluster " << i+1 << " with weight: "<< weight << std::endl;
#endif
      return 1;
    }
    return 0;
  }





};




// Allows to return class objects to R
RCPP_EXPOSED_CLASS(UserStream)
  RCPP_EXPOSED_CLASS(CF)

  // Expose class members and methods to R
  RCPP_MODULE(MOD_userStream){
    using namespace Rcpp;

    class_<UserStream>("UserStream")
      .constructor<double, double, int>()

    .field("r", &UserStream::r)
    .field("lambda", &UserStream::lambda)
    .field("tnow", &UserStream::tnow)
    .field("mean", &UserStream::mean)
    .field("dev", &UserStream::dev)
    .field("count", &UserStream::count)
    .field("currentSD", &UserStream::count)
    .field("currentMean", &UserStream::count)
    .field("numScales", &UserStream::numScales)
    .field("tgap", &UserStream::tgap)
    .field("omega", &UserStream::omega)
    .field("n", &UserStream::n)
    .field("d", &UserStream::d)

    .method("update", &UserStream::update)
    .method("get_microclusters", &UserStream::get_microclusters)
    .method("get_microweights", &UserStream::get_microweights)
    .method("getAssignment", &UserStream::getAssignment)
    .method("getMean", &UserStream::getMean)
    .method("getVariance", &UserStream::getVariance)
    .method("getSD", &UserStream::getSD)
    .method("removeInsufficient", &UserStream::removeInsufficient)
    .method("getClosest", &UserStream::getClosest)
    ;

    class_<CF>("CF")
      .constructor<Rcpp::NumericVector, Rcpp::NumericVector, double, int>()

    .field("LS", &CF::LS)
    .field("SS", &CF::SS)
    .field("t", &CF::t)
    .field("weight", &CF::weight)

    .method("getCentroid", &CF::getCentroid )
    ;

  }

