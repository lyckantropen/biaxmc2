/* 
 * File:   main.cpp
 * Author: karol
 *
 * Created on 8 grudzień 2009, 17:42
 */

#include "base.h"
#include "boost.h"
#include "std.h"
#include "valarray_external.h"
#include "SimulationDB.h"

/*
 * todo: Mathematica output
 */
typedef enum { table, mathematica, maple, count } output_t ;

std::vector<std::string> test_in(std::vector<std::string> src, std::vector<std::string> dst){
    std::vector<std::string> res;
    for(std::vector<std::string>::iterator i = src.begin();i!=src.end();i++){
        if(std::find(dst.begin(),dst.end(),*i)!=dst.end())
            res.push_back(*i);
    }
    return res;
}

int index(const std::vector<std::string> & s, const std::string & find) {
    for(int i=0;i<s.size();i++)
        if(s[i]==find) return i;
    return -1;
}


void    do_count(const std::string & data_type,const std::vector<std::string> & columns,boostbase::base & db,const boostbase::tween_t_proxy & betweens,const boostbase::pair_t_proxy & wheres){
	if(data_type=="final_properties" || data_type=="properties")
		std::cout << db.get<PRE79MeanProperties>(wheres,betweens).size() << std::endl;
	if(data_type=="lattice" || data_type=="final_lattice")
		std::cout << db.get<Lattice>(wheres,betweens).size() << std::endl;
	if(data_type=="properties_evolution" || data_type=="thermalization_history")
		std::cout << db.get<PRE79StandardProperties>(wheres,betweens).size() << std::endl;
}
void    table_output(const std::string & data_type,const std::vector<std::string> & columns,boostbase::base & db,const boostbase::tween_t_proxy & betweens,const boostbase::pair_t_proxy & wheres,bool recalculate=false){

    if(data_type=="thermalization_history") {
        std::vector<PRE79StandardProperties> whatwegot= db.get<PRE79StandardProperties>(wheres,betweens);

        std::cout << "## some of the columns may imply an additional column with error values\n";
        std::cout << "# ";

        foreach(const std::string & col,columns){
            std::cout << col << " ";
        }

        std::cout << std::endl;
        foreach(const PRE79StandardProperties & prop,whatwegot){
            for(int i=0;i<prop.GetNCycles();i++){
            foreach(const std::string & column,columns){
                if(column=="energy")
                    std::cout << prop.EnergyEvolution()[i] << "\t";
                if(column=="time")
                    std::cout << i << "\t";
                if(column=="parity")
                    std::cout << prop.ParityEvolution()[i] << "\t";

                if(column=="d200corz")
                    std::cout << prop.Delta200ZCorrelationEvolution()[i] << "\t";
                if(column=="d222corz")
                    std::cout << prop.Delta222ZCorrelationEvolution()[i] << "\t";
                if(column=="d220corz")
                    std::cout << prop.Delta220ZCorrelationEvolution()[i] << "\t";
                if(column=="d200corx")
                    std::cout << prop.Delta200XCorrelationEvolution()[i] << "\t";
                if(column=="d222corx")
                    std::cout << prop.Delta222XCorrelationEvolution()[i] << "\t";
                if(column=="d220corx")
                    std::cout << prop.Delta220XCorrelationEvolution()[i] << "\t";
                if(column=="d200cory")
                    std::cout << prop.Delta200YCorrelationEvolution()[i] << "\t";
                if(column=="d222cory")
                    std::cout << prop.Delta222YCorrelationEvolution()[i] << "\t";
                if(column=="d220cory")
                    std::cout << prop.Delta220YCorrelationEvolution()[i] << "\t";

                if(column=="d322cor")
                    std::cout << prop.Delta322CorrelationEvolution()[i] << "\t";

                if(column=="paritycor")
                    std::cout << prop.ParityCorrelationEvolution()[i] << "\t";
                }
            	std::cout << std::endl;
            }
        }
    }

    if(data_type=="final_properties" || data_type=="properties") {
        /**** wczytywanie metadanych ***/
        std::vector<std::string> m_columns = test_in(columns,db.columns());
        std::vector<std::vector<std::string> > metadata;
        if(m_columns.size())
                metadata = db.get_metadata(m_columns,wheres,betweens);
        
        std::vector<PRE79MeanProperties> whatwegot;
	if(!recalculate)
		whatwegot = db.get<PRE79MeanProperties>(wheres,betweens);
	else {
        	boostbase::pair_t_proxy nwheres = wheres;
	        nwheres.pop();
		nwheres(std::string("data_type"),std::string("properties_evolution"));
                
                std::vector<std::string> names;
                names.push_back("temperature");
                names.push_back("lambda");
                names.push_back("tau");
                names.push_back("field");
                
                std::vector<std::string> h_columns = test_in(names,db.columns());
                std::vector<std::vector<std::string> > h_data = db.get_metadata(h_columns,nwheres,betweens);
                
		//std::vector<PRE79MeanProperties> wwg = db.get<PRE79MeanProperties>(wheres,betweens);
		
                std::vector<PRE79StandardProperties> wwg2 = db.get<PRE79StandardProperties>(nwheres,betweens);
		//if(wwg2.size()<wwg.size()) std::cout << "#WARINING: insufficient data to recalculate. output may not be very satisfying\n";
                int size = wwg2.size();
                whatwegot.resize(size);

                #pragma omp parallel for ordered
		for(int i=0;i<size;i++){
                    
                        double t = std::atof(h_data[i][index(h_columns,"temperature")].c_str());
                        double l = std::atof(h_data[i][index(h_columns,"lambda")].c_str());
                        double tau = std::atof(h_data[i][index(h_columns,"tau")].c_str());
                        double h = std::atof(h_data[i][index(h_columns,"field")].c_str());
                        
                        std::cout << "t: " << t << ", l: " << ", tau: " << tau << ", h:" << h << std::endl;
                        
			//PRE79StandardHamiltonian H(wwg[i].Temperature(),wwg[i].Lambda(),wwg[i].Tau(),wwg[i].Field());
			
                        PRE79StandardHamiltonian H(t,l,tau,h);
			whatwegot[i]=PRE79MeanProperties(wwg2[i],H);
		}
	}

       
        /*	
        foreach(std::vector<std::string> & m, metadata){
            foreach(std::string & a, m){
                std::cout << a << ",";
            }
            std::cout << std::endl;
        }
	*/
        
        std::cout << "## some of the columns may imply an additional column with error values\n";
        std::cout << "# ";
        foreach(const std::string & col,columns){
            std::cout << col << " ";
        }
        std::cout << std::endl;

        int count = 0;
        foreach(PRE79MeanProperties & prop,whatwegot){
            //std::cout << prop << std::endl;
            prop.CalculateMeanTensors();
            foreach(const std::string & column,columns){
                //if(column=="temperature")
                //    std::cout << prop.Temperature() << "\t";
                if(column=="specific_heat")
                    std::cout << prop.SpecificHeat().TableForm() << "\t";
                if(column=="energy")
                    std::cout << prop.TemporalMeanEnergyPerMolecule().TableForm() << "\t";
                if(column=="parity")
                    std::cout << prop.TemporalMeanParity().TableForm() << "\t";
                if(column=="d200z_from_correlation")
                    std::cout << prop.Delta200ZByCorrelation().TableForm() << "\t";
                if(column=="d222z_from_correlation")
                    std::cout << prop.Delta222ZByCorrelation().TableForm() << "\t";
                if(column=="d200x_from_correlation")
                    std::cout << prop.Delta200XByCorrelation().TableForm() << "\t";
                if(column=="d222x_from_correlation")
                    std::cout << prop.Delta222XByCorrelation().TableForm() << "\t";
                if(column=="d200y_from_correlation")
                    std::cout << prop.Delta200YByCorrelation().TableForm() << "\t";
                if(column=="d222y_from_correlation")
                    std::cout << prop.Delta222YByCorrelation().TableForm() << "\t";
                if(column=="d322_from_correlation")
                    std::cout << prop.Delta322ByCorrelation().TableForm() << "\t";
                if(column=="parity_from_correlation")
                    std::cout << prop.ParityByCorrelation().TableForm() << "\t";
                if(column=="parity_sus")
                    std::cout << prop.ParitySusceptibility().TableForm() << "\t";
                if(column=="d200z_from_correlation_sus")
                    std::cout << prop.Delta200ZByCorrelationSusceptibility().TableForm() << "\t";
                if(column=="d222z_from_correlation_sus")
                    std::cout << prop.Delta222ZByCorrelationSusceptibility().TableForm() << "\t";
                if(column=="d200x_from_correlation_sus")
                    std::cout << prop.Delta200XByCorrelationSusceptibility().TableForm() << "\t";
                if(column=="d222x_from_correlation_sus")
                    std::cout << prop.Delta222XByCorrelationSusceptibility().TableForm() << "\t";
                if(column=="d200y_from_correlation_sus")
                    std::cout << prop.Delta200YByCorrelationSusceptibility().TableForm() << "\t";
                if(column=="d222y_from_correlation_sus")
                    std::cout << prop.Delta222YByCorrelationSusceptibility().TableForm() << "\t";
                if(column=="d322_from_correlation_sus")
                    std::cout << prop.Delta322ByCorrelationSusceptibility().TableForm() << "\t";
                if(column=="parity_from_correlation_sus")
                    std::cout << prop.ParityByCorrelationSusceptibility().TableForm() << "\t";
                /*if(column=="tau")
                    std::cout << prop.Tau() << "\t";
                if(column=="lambda")
                    std::cout << prop.Lambda() << "\t";
                //to pole nie może nazywać się h ani H, ponieważ zachodzi kolizja z wysokością siatki H
                if(column=="field")
                    std::cout << prop.Field() << "\t";
                 */
                if(column=="mean_d200corz")
                    std::cout << prop.Delta200ZMeanCorrelation() << "\t";
                if(column=="mean_d222corz")
                    std::cout << prop.Delta222ZMeanCorrelation() << "\t";
                if(column=="mean_d220corz")
                    std::cout << prop.Delta220ZMeanCorrelation() << "\t";
                if(column=="mean_d200corx")
                    std::cout << prop.Delta200XMeanCorrelation() << "\t";
                if(column=="mean_d222corx")
                    std::cout << prop.Delta222XMeanCorrelation() << "\t";
                if(column=="mean_d220corx")
                    std::cout << prop.Delta220XMeanCorrelation() << "\t";
                if(column=="mean_d200cory")
                    std::cout << prop.Delta200YMeanCorrelation() << "\t";
                if(column=="mean_d222cory")
                    std::cout << prop.Delta222YMeanCorrelation() << "\t";
                if(column=="mean_d220cory")
                    std::cout << prop.Delta220YMeanCorrelation() << "\t";
                if(column=="mean_d200")
                    std::cout << prop.MeanDelta200() << "\t";
                if(column=="mean_d220")
                    std::cout << prop.MeanDelta220() << "\t";
                if(column=="mean_d202")
                    std::cout << prop.MeanDelta202() << "\t";
                if(column=="mean_d222")
                    std::cout << prop.MeanDelta222() << "\t";


                if(column=="mean_qx")
                    std::cout << prop.MeanQxTensor() << "\t" ;
                if(column=="mean_qy")
                    std::cout << prop.MeanQyTensor() << "\t" ;
                if(column=="mean_qz")
                    std::cout << prop.MeanQzTensor() << "\t" ;

                if(column=="T20T20z")
                    std::cout << prop.MeanT20T20Z().TableForm() << "\t";
                if(column=="T20T20z_sus")
                    std::cout << prop.MeanT20T20ZSusceptibility().TableForm() << "\t";
                if(column=="T22T22z")
                    std::cout << prop.MeanT22T22Z().TableForm() << "\t";
                if(column=="T22T22z_sus")
                    std::cout << prop.MeanT22T22ZSusceptibility().TableForm() << "\t";
                if(column=="T20T20x")
                    std::cout << prop.MeanT20T20X().TableForm() << "\t";
                if(column=="T20T20x_sus")
                    std::cout << prop.MeanT20T20XSusceptibility().TableForm() << "\t";
                if(column=="T22T22x")
                    std::cout << prop.MeanT22T22X().TableForm() << "\t";
                if(column=="T22T22x_sus")
                    std::cout << prop.MeanT22T22XSusceptibility().TableForm() << "\t";
                if(column=="T20T20y")
                    std::cout << prop.MeanT20T20Y().TableForm() << "\t";
                if(column=="T20T20y_sus")
                    std::cout << prop.MeanT20T20YSusceptibility().TableForm() << "\t";
                if(column=="T22T22y")
                    std::cout << prop.MeanT22T22Y().TableForm() << "\t";
                if(column=="T22T22y_sus")
                    std::cout << prop.MeanT22T22YSusceptibility().TableForm() << "\t";

                if(column=="tetrahedral_correlation" || column=="mean_d322cor")
                    std::cout << prop.Delta322MeanCorrelation() << "\t";
                if(column=="parity_correlation")
                    std::cout << prop.ParityMeanCorrelation() << "\t";
                
                //metadane
                if(index(m_columns,column)!=-1){
                    std::cout << metadata[count][index(m_columns,column)] << "\t";
                }

            }
            std::cout << std::endl;
            count++;
        }
        
    }
    if(data_type=="lattice" || data_type=="final_lattice") {
        std::vector<Lattice> whatwegot= db.get<Lattice>(wheres,betweens);

        foreach(const Lattice & lat,whatwegot){
            std::cout << "## some of the values (e.g. vectors) may span over many columns\n";
            std::cout << "# ";
            foreach(const std::string & col,columns){
                std::cout << col << " ";
            }
            std::cout << std::endl;
            for(int p=0;p<lat.GetN();p++){
                foreach(const std::string & column,columns){
                    const Particle & cp = lat.GetParticles()[p];
                    if(column=="index")
                        std::cout << p << "\t";
                    if(column=="parity")
                        std::cout << cp.GetParity() << "\t";
                    if(column=="orientation")
                        std::cout << cp.GetX() << "\t";
                    if(column=="T")
                        std::cout << cp.GetT() << "\t";
                    if(column=="energy")
                        std::cout << cp.GetEnergy() << "\t";
                    if(column=="Ex")
                        std::cout << cp.GetEX() << "\t";
                    if(column=="Ey")
                        std::cout << cp.GetEY() << "\t";
                    if(column=="Ez")
                        std::cout << cp.GetEZ() << "\t";
                    if(column=="Qx")
                        std::cout << cp.GetQX() << "\t";
                    if(column=="Qy")
                        std::cout << cp.GetQY() << "\t";
                    if(column=="Qz")
                        std::cout << cp.GetQZ() << "\t";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
    if(data_type=="properties_evolution"){
        //TODO: tutaj jest segfault
        std::vector<PRE79StandardProperties> whatwegot= db.get<PRE79StandardProperties>(wheres,betweens);

        for(int p=0;p<whatwegot.size();p++){
            const PRE79StandardProperties & prop = whatwegot[p];
            std::cout << "## some of the values (e.g. vectors) may span over many columns\n";
            std::cout << "# ";
            foreach(const std::string & col,columns){
                std::cout << col << " ";
            }
            std::cout << std::endl;
            for(int t=0;t<prop.GetNCycles();t++){
                foreach(const std::string & column,columns){
                    if(column=="time")
                        std::cout << t << "\t";
                    if(column=="energy")
                        std::cout << prop.EnergyEvolution()[t] << "\t";
                    if(column=="parity")
                        std::cout << prop.ParityEvolution()[t] << "\t";

                    if(column=="d200corz")
                        std::cout << prop.Delta200ZCorrelationEvolution()[t] << "\t";
                    if(column=="d222corz")
                        std::cout << prop.Delta222ZCorrelationEvolution()[t] << "\t";
                    if(column=="d220corz")
                        std::cout << prop.Delta220ZCorrelationEvolution()[t] << "\t";
                    if(column=="d200corx")
                        std::cout << prop.Delta200XCorrelationEvolution()[t] << "\t";
                    if(column=="d222corx")
                        std::cout << prop.Delta222XCorrelationEvolution()[t] << "\t";
                    if(column=="d220corx")
                        std::cout << prop.Delta220XCorrelationEvolution()[t] << "\t";
                    if(column=="d200cory")
                        std::cout << prop.Delta200YCorrelationEvolution()[t] << "\t";
                    if(column=="d222cory")
                        std::cout << prop.Delta222YCorrelationEvolution()[t] << "\t";
                    if(column=="d220cory")
                        std::cout << prop.Delta220YCorrelationEvolution()[t] << "\t";

                    if(column=="d322cor")
                        std::cout << prop.Delta322CorrelationEvolution()[t] << "\t";

                    if(column=="paritycor")
                        std::cout << prop.ParityCorrelationEvolution()[t] << "\t";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }

}

/**
 * tutaj ignorujemy columns, tylko wypisujemy wszystko
 */
void mathematica_output(const std::string & data_type,const std::vector<std::string> & /*columns*/,boostbase::base & db,const boostbase::tween_t_proxy & betweens,const boostbase::pair_t_proxy & wheres,bool recalculate=false){
    if(data_type=="properties" || data_type=="final_properties"){
        std::vector<PRE79MeanProperties> whatwegot;
	if(!recalculate)
		whatwegot = db.get<PRE79MeanProperties>(wheres,betweens);
	else {
        	boostbase::pair_t_proxy nwheres = wheres;
	        nwheres.pop();
		nwheres(std::string("data_type"),std::string("properties_evolution"));
		std::vector<PRE79MeanProperties> wwg = db.get<PRE79MeanProperties>(wheres,betweens);
		std::vector<PRE79StandardProperties> wwg2 = db.get<PRE79StandardProperties>(nwheres,betweens);
		if(wwg2.size()<wwg.size()) std::cout << "#WARINING: insufficient data to recalculate. output may not be very satisfying\n";
                int size = wwg2.size();
                whatwegot.resize(size);

                #pragma omp parallel for ordered
		for(int i=0;i<size;i++){
			PRE79StandardHamiltonian H(wwg[i].Temperature(),wwg[i].Lambda(),wwg[i].Tau(),wwg[i].Field());
			whatwegot[i]=PRE79MeanProperties(wwg2[i],H);
		}
	}
        vect temperature(0.0,whatwegot.size());
        //foreach(PRE79MeanProperties & prop,whatwegot){
        for(int i=0;i<whatwegot.size();i++){
            PRE79MeanProperties & prop = whatwegot[i];
            prop.CalculateMeanTensors();
            //współrzędne - wszystko jest w zależności od T,h,lambda i tau
            std::stringstream coord;
            temperature[i]=prop.Temperature();
            coord << "[" << prop.Temperature() << "," << prop.Field() << "," << prop.Lambda() << "," << prop.Tau() << "]" ;

            std::cout << "SpecificHeat" << coord.str() << "=" << prop.SpecificHeat().MathematicaForm() << ";\n";
            std::cout << "Fluctuation" << coord.str() << "=" << prop.Fluctuation().MathematicaForm() << ";\n";
            std::cout << "TemporalMeanEnergyPerMolecule" << coord.str() << "=" << prop.TemporalMeanEnergyPerMolecule().MathematicaForm() << ";\n";
            std::cout << "TemporalMeanParity" << coord.str() << "=" << prop.TemporalMeanParity().MathematicaForm() << ";\n";
            std::cout << "ParitySusceptibility" << coord.str() << "=" << prop.ParitySusceptibility().MathematicaForm() << ";\n";
            std::cout << "Delta200ZByCorrelation" << coord.str() << "=" << prop.Delta200ZByCorrelation().MathematicaForm() << ";\n";
            std::cout << "Delta222ZByCorrelation" << coord.str() << "=" << prop.Delta222ZByCorrelation().MathematicaForm() << ";\n";
            std::cout << "Delta200XByCorrelation" << coord.str() << "=" << prop.Delta200XByCorrelation().MathematicaForm() << ";\n";
            std::cout << "Delta222XByCorrelation" << coord.str() << "=" << prop.Delta222XByCorrelation().MathematicaForm() << ";\n";
            std::cout << "Delta200YByCorrelation" << coord.str() << "=" << prop.Delta200YByCorrelation().MathematicaForm() << ";\n";
            std::cout << "Delta222YByCorrelation" << coord.str() << "=" << prop.Delta222YByCorrelation().MathematicaForm() << ";\n";

            std::cout << "Delta322ByCorrelation" << coord.str() << "=" << prop.Delta322ByCorrelation().MathematicaForm() << ";\n";
            std::cout << "ParityByCorrelation" << coord.str() << "=" << prop.ParityByCorrelation().MathematicaForm() << ";\n";

            std::cout << "Delta200ZByCorrelationSusceptibility" << coord.str() << "=" << prop.Delta200ZByCorrelationSusceptibility().MathematicaForm() << ";\n";
            std::cout << "Delta222ZByCorrelationSusceptibility" << coord.str() << "=" << prop.Delta222ZByCorrelationSusceptibility().MathematicaForm() << ";\n";
            std::cout << "Delta200XByCorrelationSusceptibility" << coord.str() << "=" << prop.Delta200XByCorrelationSusceptibility().MathematicaForm() << ";\n";
            std::cout << "Delta222XByCorrelationSusceptibility" << coord.str() << "=" << prop.Delta222XByCorrelationSusceptibility().MathematicaForm() << ";\n";
            std::cout << "Delta200YByCorrelationSusceptibility" << coord.str() << "=" << prop.Delta200YByCorrelationSusceptibility().MathematicaForm() << ";\n";
            std::cout << "Delta222YByCorrelationSusceptibility" << coord.str() << "=" << prop.Delta222YByCorrelationSusceptibility().MathematicaForm() << ";\n";

            std::cout << "Delta322ByCorrelationSusceptibility" << coord.str() << "=" << prop.Delta322ByCorrelationSusceptibility().MathematicaForm() << ";\n";
            std::cout << "ParityByCorrelationSusceptibility" << coord.str() << "=" << prop.ParityByCorrelationSusceptibility().MathematicaForm() << ";\n";

            std::cout << "MeanT20T20Z" << coord.str() << "=" << prop.MeanT20T20Z().MathematicaForm() << ";\n";
            std::cout << "MeanT20T20ZSusceptibility" << coord.str() << "=" << prop.MeanT20T20ZSusceptibility().MathematicaForm() << ";\n";
            std::cout << "MeanT22T22Z" << coord.str() << "=" << prop.MeanT22T22Z().MathematicaForm() << ";\n";
            std::cout << "MeanT22T22ZSusceptibility" << coord.str() << "=" << prop.MeanT22T22ZSusceptibility().MathematicaForm() << ";\n";

            std::cout << "MeanT20T20X" << coord.str() << "=" << prop.MeanT20T20X().MathematicaForm() << ";\n";
            std::cout << "MeanT20T20XSusceptibility" << coord.str() << "=" << prop.MeanT20T20XSusceptibility().MathematicaForm() << ";\n";
            std::cout << "MeanT22T22X" << coord.str() << "=" << prop.MeanT22T22X().MathematicaForm() << ";\n";
            std::cout << "MeanT22T22XSusceptibility" << coord.str() << "=" << prop.MeanT22T22XSusceptibility().MathematicaForm() << ";\n";

            std::cout << "MeanT20T20Y" << coord.str() << "=" << prop.MeanT20T20Y().MathematicaForm() << ";\n";
            std::cout << "MeanT20T20YSusceptibility" << coord.str() << "=" << prop.MeanT20T20YSusceptibility().MathematicaForm() << ";\n";
            std::cout << "MeanT22T22Y" << coord.str() << "=" << prop.MeanT22T22Y().MathematicaForm() << ";\n";
            std::cout << "MeanT22T22YSusceptibility" << coord.str() << "=" << prop.MeanT22T22YSusceptibility().MathematicaForm() << ";\n";


            std::cout << "Delta200ZMeanCorrelation" << coord.str() << "=" << MathematicaForm(prop.Delta200ZMeanCorrelation()) << ";\n";
            std::cout << "Delta220ZMeanCorrelation" << coord.str() << "=" << MathematicaForm(prop.Delta220ZMeanCorrelation()) << ";\n";
            std::cout << "Delta222ZMeanCorrelation" << coord.str() << "=" << MathematicaForm(prop.Delta222ZMeanCorrelation()) << ";\n";

            std::cout << "Delta200XMeanCorrelation" << coord.str() << "=" << MathematicaForm(prop.Delta200XMeanCorrelation()) << ";\n";
            std::cout << "Delta220XMeanCorrelation" << coord.str() << "=" << MathematicaForm(prop.Delta220XMeanCorrelation()) << ";\n";
            std::cout << "Delta222XMeanCorrelation" << coord.str() << "=" << MathematicaForm(prop.Delta222XMeanCorrelation()) << ";\n";

            std::cout << "Delta200YMeanCorrelation" << coord.str() << "=" << MathematicaForm(prop.Delta200YMeanCorrelation()) << ";\n";
            std::cout << "Delta220YMeanCorrelation" << coord.str() << "=" << MathematicaForm(prop.Delta220YMeanCorrelation()) << ";\n";
            std::cout << "Delta222YMeanCorrelation" << coord.str() << "=" << MathematicaForm(prop.Delta222YMeanCorrelation()) << ";\n";

            std::cout << "Delta200" << coord.str() << "=" << prop.MeanDelta200() << ";\n";
            std::cout << "Delta220" << coord.str() << "=" << prop.MeanDelta220() << ";\n";
            std::cout << "Delta202" << coord.str() << "=" << prop.MeanDelta202() << ";\n";
            std::cout << "Delta222" << coord.str() << "=" << prop.MeanDelta222() << ";\n";

            std::cout << "MeanQxTensor" << coord.str() << "=" << MathematicaForm(prop.MeanQxTensor()) << ";\n";
            std::cout << "MeanQyTensor" << coord.str() << "=" << MathematicaForm(prop.MeanQyTensor()) << ";\n";
            std::cout << "MeanQzTensor" << coord.str() << "=" << MathematicaForm(prop.MeanQzTensor()) << ";\n";


        }
        std::cout << "Temp=" << MathematicaForm(temperature) << ";\n";
    }
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
            for(int l=0;l<lat.GetL();l++)
             for(int w=0;w<lat.GetW();w++)
              for(int h=0;h<lat.GetH();h++){
                  int p=h*(lat.GetL()*lat.GetW())+l*lat.GetW() + w;
                  std::stringstream coord;
                  coord << "[" << prop.Temperature() << "," << prop.Field() << "," << prop.Lambda() << "," << prop.Tau() ;
                  coord << "," << l << "," << w << "," << h <<"]";
                  const Particle & cp = lat.GetParticles()[p];

                  std::cout << "Parity" << coord.str() << "=" << cp.GetParity() << ";\n";
                  std::cout << "Energy" << coord.str() << "=" << cp.GetEnergy() << ";\n";
                  std::cout << "X" << coord.str() << "=" << MathematicaForm(cp.GetX()) << ";\n";
                  std::cout << "T" << coord.str() << "=" << MathematicaForm(cp.GetT()) << ";\n";
                  std::cout << "Ex" << coord.str() << "=" << MathematicaForm(cp.GetEX()) << ";\n";
                  std::cout << "Ey" << coord.str() << "=" << MathematicaForm(cp.GetEY()) << ";\n";
                  std::cout << "Ez" << coord.str() << "=" << MathematicaForm(cp.GetEZ()) << ";\n";
                  std::cout << "Qx" << coord.str() << "=" << MathematicaForm(cp.GetQX()) << ";\n";
                  std::cout << "Qy" << coord.str() << "=" << MathematicaForm(cp.GetQY()) << ";\n";
                  std::cout << "Qz" << coord.str() << "=" << MathematicaForm(cp.GetQZ()) << ";\n";
              }
        }
    }
}


int main(int argc, char** argv)
{
    std::string                 data_type="final_properties";
    std::string dbfile="";
    std::string dbdir="";
    boostbase::tween_t_proxy    betweens;
    boostbase::pair_t_proxy     wheres;
    std::vector<std::string>    columns;
    bool sqlite_debug=false;
    bool recalculate=false;
    output_t output_type = table;

    for(int i=1;i<argc;i++){
        if(std::string(argv[i])=="--db"){
            dbfile=argv[i+1];
            dbdir=argv[i+2];
        }
        if(std::string(argv[i])=="--data-type")
            data_type=argv[i+1];
        if(std::string(argv[i])=="--columns"){
            std::string columns_string = argv[i+1];
            boost::split(columns,columns_string,boost::is_any_of(","));
        }
        if(std::string(argv[i])=="--between")
            betweens(std::string(argv[i+1]),std::atof(argv[i+2]),std::atof(argv[i+3]));
        if(std::string(argv[i])=="--where")
            wheres(std::string(argv[i+1]),std::string(argv[i+2]));
        if(std::string(argv[i])=="--sqlite-debug")
            sqlite_debug=true;
        if(std::string(argv[i])=="--day"){
            pt::ptime day_start = pt::time_from_string(std::string(argv[i+1])+std::string(" 00:00:00.000000"));
            pt::ptime day_end = day_start + pt::hours(24);
            betweens(std::string("date"),pt::to_simple_string(day_start),pt::to_simple_string(day_end));
        }
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
        if(std::string(argv[i])=="--yesterday"){
            pt::ptime day_end = pt::second_clock::local_time();
            day_end-=day_end.time_of_day();
            pt::ptime day_start = day_end - pt::hours(24);
            betweens(std::string("date"),pt::to_simple_string(day_start),pt::to_simple_string(day_end));
        }
        if(std::string(argv[i])=="--last-week"){
            pt::ptime day_end = pt::second_clock::local_time();
            pt::ptime day_start = day_end - pt::hours(7*24);
            betweens(std::string("date"),pt::to_simple_string(day_start),pt::to_simple_string(day_end));
        }

        if(std::string(argv[i])=="--days-ago"){
            pt::ptime day_end = pt::second_clock::local_time();
            pt::ptime day_start = day_end - pt::hours(std::atof(argv[i+1])*24);
            betweens(std::string("date"),pt::to_simple_string(day_start),pt::to_simple_string(day_end));
        }
        if(std::string(argv[i])=="--hours-ago"){
            pt::ptime day_end = pt::second_clock::local_time();
            pt::ptime day_start = day_end - pt::hours(std::atof(argv[i+1]));
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
        if(std::string(argv[i])=="--mathematica")
            output_type=mathematica;
	if(std::string(argv[i])=="--count")
	    output_type=count;
	if(std::string(argv[i])=="--recalculate")
	    recalculate=true;

    }
    if(dbfile=="" || dbdir==""){
        std::cout << "please specify database file and directory with --db\n";
        std::exit(1);
    }

    wheres(std::string("data_type"),data_type);
    std::cout << std::setprecision(12) << std::fixed ;

    //readonly 
    boostbase::base db(dbfile,dbdir,true);
   

    switch(output_type){
        case table:
            table_output(data_type,columns,db,betweens,wheres,recalculate);
            break;
        case mathematica:
            mathematica_output(data_type,columns,db,betweens,wheres,recalculate);
            break;
	case count:
	    do_count(data_type,columns,db,betweens,wheres);
	    break;
    }

    if(sqlite_debug)
        std::cout << std::endl << db.log().str() << std::endl;

    return (EXIT_SUCCESS);
}

