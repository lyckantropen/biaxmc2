/* 
 * File:   main.cpp
 * Author: karol
 *
 * Created on 8 grudzień 2009, 17:42
 */

#include "base.h"
#include "boost.h"
#include "std.h"
#include "SimulationDB.h"


int main(int argc, char** argv)
{
    std::string                 data_type="final_lattice";
    std::string dbfile="";
    std::string dbdir="";
    boostbase::tween_t_proxy    betweens;
    boostbase::pair_t_proxy     wheres;
    std::vector<std::string>    columns;
    bool sqlite_debug=false;

    for(int i=1;i<argc;i++){
        if(std::string(argv[i])=="--db"){
            dbfile=argv[i+1];
            dbdir=argv[i+2];
        }
        if(std::string(argv[i])=="--day"){
            pt::ptime day_start = pt::time_from_string(std::string(argv[i+1])+std::string(" 00:00:00.000000"));
            pt::ptime day_end = day_start + pt::hours(24);
            betweens(std::string("date"),pt::to_simple_string(day_start),pt::to_simple_string(day_end));
        }
        if(std::string(argv[i])=="--sqlite-debug")
            sqlite_debug=true;
        if(std::string(argv[i])=="--today"){
            pt::ptime day_start = pt::second_clock::local_time();
            day_start-=day_start.time_of_day();
            pt::ptime day_end = day_start + pt::hours(24);
            betweens(std::string("date"),pt::to_simple_string(day_start),pt::to_simple_string(day_end));
        }
        if(std::string(argv[i])=="--days"){
            pt::ptime day_start = pt::time_from_string(std::string(argv[i+1])+std::string(" 00:00:00.000000"));
            pt::ptime day_end = pt::time_from_string(std::string(argv[i+2])+std::string(" 00:00:00.000000"));
            betweens(std::string("date"),pt::to_simple_string(day_start),pt::to_simple_string(day_end));
        }
        if(std::string(argv[i])=="--between")
            betweens(std::string(argv[i+1]),std::atof(argv[i+2]),std::atof(argv[i+3]));
        if(std::string(argv[i])=="--where")
            wheres(std::string(argv[i+1]),std::string(argv[i+2]));
        if(std::string(argv[i])=="--yesterday"){
            pt::ptime day_end = pt::second_clock::local_time();
            day_end-=day_end.time_of_day();
            pt::ptime day_start = day_end - pt::hours(24);
            betweens(std::string("date"),pt::to_simple_string(day_start),pt::to_simple_string(day_end));
        }
        if(std::string(argv[i])=="--hours"){
            pt::ptime hour_start;
            pt::ptime hour_end;

            pt::ptime day_start = pt::second_clock::local_time();
            day_start-=day_start.time_of_day();

            hour_start = day_start + pt::duration_from_string(std::string(argv[i+1]));
            hour_end = day_start + pt::duration_from_string(std::string(argv[i+2]));

            betweens(std::string("date"),pt::to_simple_string(hour_start),pt::to_simple_string(hour_end));
        }

    }
    if(dbfile=="" || dbdir==""){
        std::cout << "please specify database file and directory with --db\n";
        std::exit(1);
    }

    wheres(std::string("data_type"),data_type);
    std::cout << std::setprecision(12) << std::fixed ;

    //readonly 
    boostbase::base db(dbfile,dbdir,true);

    if(data_type=="lattice" || data_type=="final_lattice") {
        std::vector<Lattice> whatwegot= db.get<Lattice>(wheres,betweens);
        boostbase::pair_t_proxy nwheres = wheres;
        nwheres.pop();
        nwheres(std::string("data_type"),std::string("final_properties"));
        std::vector<PRE79MeanProperties> props= db.get<PRE79MeanProperties>(nwheres,betweens);

        //foreach(const Lattice & lat,whatwegot) {
        for(int i=0;i<whatwegot.size();i++) {
            const Lattice & lat = whatwegot[i];
            const PRE79MeanProperties & prop = props[i];

            std::stringstream coord;
            coord << "[" << prop.Temperature() << "," << prop.Field() << "," << prop.Lambda() << "," << prop.Tau() << "]" ;

	    vect MeanQxTensor(0.0,6);
	    vect MeanQyTensor(0.0,6);
	    vect MeanQzTensor(0.0,6);
	    vect MeanQ2(0.0,6);
	    vect MeanQ3(0.0,6);

            for(int l=0;l<lat.GetL();l++)
             for(int w=0;w<lat.GetW();w++)
              for(int h=0;h<lat.GetH();h++){
                  int p=h*(lat.GetL()*lat.GetW())+l*lat.GetW() + w;
                  const Particle & cp = lat.GetParticles()[p];

		  MeanQxTensor+=cp.GetQX();
		  MeanQyTensor+=cp.GetQY();
		  MeanQzTensor+=cp.GetQZ();

		  //vect Q1(0.0,6);
		  //vect Q2(0.0,6);
		  vect Q(0.0,6);
		  vect Q2(0.0,6);
		  Q = std::sqrt(1.5)*cp.GetQZ()-1./std::sqrt(6.)*Identity(3)+prop.Lambda()*(cp.GetQX()-cp.GetQY());
		  //Q2 = std::sqrt(1.5)*cp.GetQY()-1./std::sqrt(6.)*Identity(3)+prop.Lambda()*(cp.GetQZ()-cp.GetQX());
		  //Q1 = std::sqrt(1.5)*cp.GetQX()-1./std::sqrt(6.)*Identity(3)+prop.Lambda()*(cp.GetQY()-cp.GetQZ());

		  Q2[0]=Q[0]*Q[0] + Q[1]*Q[1] + Q[2]*Q[2];
		  Q2[1]=Q[0]*Q[1] + Q[1]*Q[3] + Q[2]*Q[4];
		  Q2[2]=Q[0]*Q[2] + Q[1]*Q[4] + Q[2]*Q[5];
		  Q2[3]=Q[1]*Q[1] + Q[3]*Q[3] + Q[4]*Q[4];
		  Q2[4]=Q[1]*Q[2] + Q[3]*Q[4] + Q[4]*Q[5];
		  Q2[5]=Q[2]*Q[2] + Q[4]*Q[4] + Q[5]*Q[5];

		  MeanQ2+=Q2;

		  MeanQ3[0]+=Q[0]*Q2[0] + Q[1]*Q2[1] + Q[2]*Q2[2];
		  MeanQ3[1]+=Q[0]*Q2[1] + Q[1]*Q2[3] + Q[2]*Q2[4];	
		  MeanQ3[2]+=Q[0]*Q2[2] + Q[1]*Q2[4] + Q[2]*Q2[5];
		  MeanQ3[3]+=Q[1]*Q2[1] + Q[3]*Q2[3] + Q[4]*Q2[4];
		  MeanQ3[4]+=Q[1]*Q2[2] + Q[3]*Q2[4] + Q[4]*Q2[5];
		  MeanQ3[5]+=Q[2]*Q2[2] + Q[4]*Q2[4] + Q[5]*Q2[5];
/*
		  MeanTrQ2[0]+=Q1[0]*Q1[0] + 2.*Q1[1]*Q1[1] + 2.*Q1[2]*Q1[2] + Q1[3]*Q1[3] + 2.*Q1[4]*Q1[4] + Q1[5]*Q1[5];
		  MeanTrQ2[1]+=Q2[0]*Q2[0] + 2.*Q2[1]*Q2[1] + 2.*Q2[2]*Q2[2] + Q2[3]*Q2[3] + 2.*Q2[4]*Q2[4] + Q2[5]*Q2[5];
		  MeanTrQ2[2]+=Q3[0]*Q3[0] + 2.*Q3[1]*Q3[1] + 2.*Q3[2]*Q3[2] + Q3[3]*Q3[3] + 2.*Q3[4]*Q3[4] + Q3[5]*Q3[5];

		  MeanTrQ3[0]+=Q1[0]*Q1[0]*Q1[0] + 3.*Q1[0]*(Q1[1]*Q1[1]+Q1[2]*Q1[2]) + Q1[3]*Q1[3]*Q1[3] + 3.*Q1[1]*Q1[1]*Q1[3] + 6.*Q1[1]*Q1[2]*Q1[4] + 3.*Q1[3]*Q1[4]*Q1[4] +
			  3.*(Q1[2]*Q1[2]+Q1[4]*Q1[4])*Q1[5] + Q1[5]*Q1[5]*Q1[5];

		  MeanTrQ3[1]+=Q2[0]*Q2[0]*Q2[0] + 3.*Q2[0]*(Q2[1]*Q2[1]+Q2[2]*Q2[2]) + Q2[3]*Q2[3]*Q2[3] + 3.*Q2[1]*Q2[1]*Q2[3] + 6.*Q2[1]*Q2[2]*Q2[4] + 3.*Q2[3]*Q2[4]*Q2[4] +
			  3.*(Q2[2]*Q2[2]+Q2[4]*Q2[4])*Q2[5] + Q2[5]*Q2[5]*Q2[5];

		  MeanTrQ3[2]+=Q3[0]*Q3[0]*Q3[0] + 3.*Q3[0]*(Q3[1]*Q3[1]+Q3[2]*Q3[2]) + Q3[3]*Q3[3]*Q3[3] + 3.*Q3[1]*Q3[1]*Q3[3] + 6.*Q3[1]*Q3[2]*Q3[4] + 3.*Q3[3]*Q3[4]*Q3[4] +
			  3.*(Q3[2]*Q3[2]+Q3[4]*Q3[4])*Q3[5] + Q3[5]*Q3[5]*Q3[5];
*/

              }
	    MeanQxTensor/=lat.GetN();
	    MeanQyTensor/=lat.GetN();
	    MeanQzTensor/=lat.GetN();

	    MeanQ2/=lat.GetN();
	    MeanQ3/=lat.GetN();

	    std::cout << "MeanQxTensor" << coord.str() << "=" << MathematicaForm(MeanQxTensor) << ";\n";
	    std::cout << "MeanQyTensor" << coord.str() << "=" << MathematicaForm(MeanQyTensor) << ";\n";
	    std::cout << "MeanQzTensor" << coord.str() << "=" << MathematicaForm(MeanQzTensor) << ";\n";
	    std::cout << "TrMeanQ2" << coord.str() << "=" << (MeanQ2[0]+MeanQ2[3]+MeanQ2[5]) << ";\n";
	    std::cout << "TrMeanQ3" << coord.str() << "=" << (MeanQ3[0]+MeanQ3[3]+MeanQ3[5]) << ";\n";

        }
    }

    if(sqlite_debug)
        std::cout << std::endl << db.log().str() << std::endl;

    return (EXIT_SUCCESS);
}

