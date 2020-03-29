package PoolHap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

//import org.apache.spark.ml.regression.LinearRegression;
//import org.apache.spark.ml.regression.LinearRegressionModel;
//import org.apache.spark.ml.regression.LinearRegressionTrainingSummary;
//import org.apache.spark.mllib.linalg.Vectors;
//import org.apache.spark.sql.Dataset;
//import org.apache.spark.sql.Row;
//import org.apache.spark.sql.SparkSession;

public class HapLASSO {
    /**
     *  @author Quan Long, 2019-03-12.
     *
     *  Given a HapConfig as potential set of input haplotypes,
     *  figure out their frequencies (with most of the being zero) using regularization
     *  The algorithm is inspired by Zhou, Zhang, Yang BIOINFORMATICS 2019 (CSHAP)
     *
     *  This can be used for the identification of either in-pool frequencies or global frequencies;
     *
     *  Regarding to LD constraints.
     *    (1) The LD constraints could be estimated statistically in a correlation framework or
     *        directed counted from reads.
     *        (1.1) in-pool LD can only be estimated by counting reads
     *              ("calcualte_LD_matrix_single_pool(String vef_file)").
     *        (1.2) global LD can be estimated statistically
     *              ("calcualte_LD_matrix_multi_pool_stat()")
     *              or by aggregating all in-pool read-counting-based LDs
     *              ("calcualte_LD_matrix_multi_pool_reads(String[] vef_file)").
     *    (2) Which LD constraints to use
     *        (2.1) One may use pool-specific LD to estimate in-pool frequencies if it is believed
     *              that within-tissue evolution matters;
     *        (2.2) Or one may use global LD to estimate in-pool frequencies if it is believed that
     *              no within-tissue evolution;
     *        (2.3) Of course, only global LD can be used to estimate global frequencies.
     *
     *  The code will write 1, H, H^H to a file in libsvm format:
     *
     *  label index1:value1 index2:value2 ...
     *  (here, in our simplified format, label is Y, the target value, index_i will be just i).
     *
     *  The first line is 1 (with the first column 1);
     *  The second block is H (with the first column MAF);
     *  The third block is H^H (with the first column an estimate of LD);
     */

    /*
     *  Input data.
     */
    HapConfig potential_haps; // all potential haplotypes out of GC (or other algorithms)

    /*
     *  Input parameters.
     */
    public double lambda; // the parameter that sets the weight of regularizer; not used in v8

    // The cutoff for a haplotype in HapConfig to be passed to the regression.
    double raw_hap_freq_cutoff;
    String memory_usage; // in bits

    // Which pool (in the HapConfig potential_haps) to work with.
    int pool_index; // WILL WORK ON GLOBAL FREQUENCES IF pool_index IS SET TO -1
    String pool_name;

    /*
     *  Calculated parameters.
     */
    double sum1_weight; // the weight for the first sample

    // The weight for each MAF sample: dimension: q, where q is the number of loci.
    double[] maf_weights;

    // The weight for each LD sample: the real dimension is: q(q - 1) / 2, where q is the number of
    // loci. For simplicity of coding, ld_weights still claims a q x q matrix, but only the off-diag
    // will be used (the rest of the entries are 0).
    double[][] ld_weights;
    int num_loci;
    double total_ld_covered_reads_counts;
    double total_ld_covered_paired_alleles;
    int[][] total_counts;
    double lasso_coverage_weight;
    double lasso_distance_max_weight;
    String sequencing_technology;
    int l0l1_regression_level;
    double lasso_cross_region_weight = 1.0;
    String species;

    /*
     *  Output and intermediate data.
     */
    // The location of sites: for fast mapping of indices when parsing the VEF files.
    HashMap<Integer, Integer> site_locations2index;
    double[] maf; // will be calculated differently for in-pool MAF and population MAF
    double[][] ld_matrix; // will be calculated differently for in-pool LD and population LD
    String prefix;
    HashMap<Integer, Integer> index_2_pos_dict;
    

    // TODO: [LEFTOVER]
    // // The indices (in Hapconfig potential_haps) of haps that enter LASSO
    // int[] quanlified_hap_index;

    // Used to assess the success/failure of LASSO. If failure, adjust lambda and try again.
    //public double r2;

    // The countcome, i.e., the coefficients from LASSO, corresponding to the haps in
    // quanlified_hap_index.
    double[] out_hap_freqs;
    String [] rscript_files;
    HashMap<Integer,String>  pos_come_from_region;
    int num_regions_from;

    /**
     *  Constructor
     *
     *  If the haplotype frequencies are not reliable (e.g. out of GC), we should set
     *  gc_hap_freq_cutoff = 0;
     *
     *  @param pool_index
     *  @param pool_name
     *  @param lambda
     *  @param potential_haps
     *  @param raw_hap_freq_cutoff
     *  @param memory_usage
     *  @param prefix  // folder path + name of the pool.
     */
    public HapLASSO(
        int pool_index,
        String pool_name,
        double lambda,
        HapConfig potential_haps,
        double raw_hap_freq_cutoff,
        String memory_usage,
        String prefix,
        double weight,
        double max_weight
        ) {
    	
    	this.lasso_coverage_weight= weight;
    	this.lasso_distance_max_weight= max_weight;
        this.pool_index = pool_index; // note: if pool_index==-1, then work on global frequencies
        this.pool_name=pool_name;
        this.lambda = lambda;
        this.raw_hap_freq_cutoff = raw_hap_freq_cutoff;
        this.potential_haps = potential_haps;
        this.num_loci = this.potential_haps.num_loci;
        // this.maf and this.ld_matrix will be created when calculating them.

        this.memory_usage = memory_usage;
        this.prefix = prefix;
       
        
        // this.potential_haps.num_global_hap and this.quanlified_hap_index will be set up in the
        // function write_H_and_1_to_file();
    }
    
    
   
    
    
    
    public HapLASSO(
    		String [] pool_names,
    		int [][] reg_region,
            String r_prefix,
            String prefix,
            int level, 
            int maxSuppSize,
            int minSize,
            double gammaMin,
            double gammaMax,
            int nGamma,
            double maf_weights, 
            String sp) throws IOException  {
    		this.species = sp;
    		int  rscript_files_len=0;
    		for (int r=0; r< reg_region.length;r++) {
    			if ((reg_region[r][1]!=-1 )  && (reg_region[r][0]!=-1 )) {
    				rscript_files_len++;
    			}
    		}
    		if (rscript_files_len==1) {
    			maxSuppSize=(int) ((double)maxSuppSize*1.0);
    		}
    		this.rscript_files= new String [rscript_files_len];
        	for (int r=0; r< reg_region.length;r++) {
        		if ((reg_region[r][1]!=-1 )  && (reg_region[r][0]!=-1 )) {
		    		for (int p=0;p< pool_names.length;p++) {
		    			String rfile = r_prefix+"/regression_level_"+ Integer.toString(level)+"_region_"+ 
		        				Integer.toString(r+1)+"_pool_" + pool_names[p]+".R";
		        		this.rscript_files[r] = rfile;
			    		BufferedWriter bw = new BufferedWriter(new FileWriter( rfile, true));
			    		bw.write("library(\'L0Learn\')\n");
			    		bw.write("tryCatch({\n");
			    		bw.write("maxSize= "+Integer.toString(maxSuppSize)+"\n");
			    		bw.write("minSize= "+Integer.toString(minSize)+"\n");
			    		bw.write("mum_sites="+Integer.toString(reg_region[r][1]-reg_region[r][0] +1)+"\n");
			    		bw.write("maf_weights= "+Double.toString(maf_weights)+"\n");
		    			String tmp ="txt<-as.matrix(read.table(\'";
		    			String infile = prefix+ pool_names[p]+"/regression_level_" + Integer.toString(level)+
		    					"_region_"+ Integer.toString(r+1)+ ".regression_in";
		    			tmp=tmp+infile+"\',sep=\' \'))";
		    			bw.write(tmp+"\n");
		    			bw.write("y<-as.numeric(txt[,2])\n");//apply(X,2,as.numeric)
		    			bw.write("X<-apply(txt[,3:ncol(txt)], 2,as.numeric)  \n");
		    			bw.write("if ((ncol(X))>minSize){\n");
	//	    			cvfit = L0Learn.cvfit(X, y, nFolds=5, seed=1, penalty="L0L1", nGamma=5, gammaMin=0.001, gammaMax=0.1, maxSuppSize=maxSize)
		    			tmp= "cvfit = L0Learn.cvfit(X, y, nFolds=5, seed=1, penalty=\"L0L1\", "
		    					+ "nGamma="+Integer.toString(nGamma)+", gammaMin="+ Double.toString(gammaMin)+","
		    					+ " gammaMax="+ Double.toString(gammaMax)+ ", maxSuppSize=maxSize"+", intercept= FALSE)";
//		    			tmp= "cvfit = L0Learn.cvfit(X, y, nFolds=5, seed=1, penalty=\"L0L1\", "
//		    					+ "nGamma="+Integer.toString(nGamma)+", gammaMin="+ Double.toString(gammaMin)+","
//		    					+ " gammaMax="+ Double.toString(gammaMax)+ ", maxSuppSize=maxSize"+", intercept= FALSE"
//		    							+ ",maxIters = 5)";
		    			bw.write(tmp+"\n");
		    			bw.write("optimalGammaIndex= which(unlist(lapply(cvfit$cvMeans, min)) == "
		    					+ "min(unlist(lapply(cvfit$cvMeans, min))))\n");
		    			bw.write("optimalLambdaIndex = which.min(cvfit$cvMeans[[optimalGammaIndex]])\n");
		    			bw.write("optimalLambda = cvfit$fit$lambda[[optimalGammaIndex]][optimalLambdaIndex]\n");
		    			if (!this.species.equals("bacteria")) {
		    				bw.write("tmp =coef(cvfit, lambda=optimalLambda, gamma=cvfit$fit$gamma[optimalGammaIndex])\n");
		    			}else {
		    				bw.write("tmp =coef(cvfit, lambda=optimalLambda, gamma=cvfit$fit$gamma[optimalGammaIndex])\n");
		    			}
		    			
		    			String outfile= prefix+ pool_names[p]+"/regression_level_" + Integer.toString(level)+
		    					"_region_"+ Integer.toString(r+1)+ ".regression_out";
		    			tmp= "write( paste( \"#Hap_ID\", \"Freq\" ,\"Haplotype\",sep = \"\\t\" ) "
		    					+ ",file=\""+outfile+"\",append=TRUE)";
		    			bw.write(tmp+"\n");
		    			bw.write("for(i in 1:length(as.vector(tmp@x))){\n");
		    			
		    			bw.write("hap = gsub(\", \",\"\",toString(round(X[, "
		    					+ "as.vector(tmp@i)[i]+1][2:(1+ mum_sites)]/maf_weights)))\n");
		    			
		    					
		    			
		    			tmp= "write(paste( paste( \"h\", toString(as.vector(tmp@i)[i]),sep = \"\" ), "
		    					+ "as.vector(tmp@x)[i],hap, sep = \"\\t\"),file=\""+outfile
		    					+"\",append=TRUE)";
		    			bw.write(tmp+"\n");
		    			bw.write("}}else {\n");
		    			tmp= "write( paste( \"#Hap_ID\", \"Freq\" ,\"Haplotype\",sep = \"\\t\" ) "
		    					+ ",file=\""+outfile+"\",append=TRUE)";
		    			bw.write(tmp+"\n");
		    			bw.write("for(i in 1:ncol(X)){\n");
		    			bw.write("hap = gsub(\", \",\"\",toString(round(X[, i][2:(1+ mum_sites)]/maf_weights)))\n");
		    			tmp= "write(paste( paste( \"h\", toString(i),sep = \"\" ), "
		    					+ "1.0/as.double(ncol(X)), hap, sep = \"\\t\"),file=\""+outfile
		    					+"\",append=TRUE)\n";
		    			bw.write(tmp);
		    			bw.write("}}\n");
		    			bw.write("print (\"FINISH\")\n");
		    			
		    			bw.write("} , error = function(err.msg){\n");
		    			bw.write("write(toString(err.msg),\""+outfile+".log\""+  ", append=TRUE)\n");
		    			bw.write("}\n");
		    			bw.write(")\n");
			    		bw.close();
		    		}
		    		
	        	}
        	}
        }
    
    public void run_L0L1learn(String [] pool_names, String rscript_path , 
    		String rscript_prefix, String outfile_prefix, 
    		int level, int region , int num_threads)  throws IOException, InterruptedException  {
    	
    	DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
    	int num_run = pool_names.length/ num_threads;
    	if( (pool_names.length% num_threads) !=0) {
    		num_run++;
    	}
    	
    	String [][] mul_pool_names = new String [ num_run][num_threads];
    	for (int n=0; n< num_run ;n++ ) {
    		for (int t=0; t< num_threads ;t++ ) {
    			if ((n*num_threads+ t  )< pool_names.length) {
    				mul_pool_names[n][t]= pool_names [n*num_threads+ t ] ;
    			}else {
    				mul_pool_names[n][t]= "NULL" ;
    			}
    		}
    	}
    	
    	for (int n=0; n< num_run ;n++ ) {
    		for (int t=0; t<num_threads;t ++ ) {
    			System.out.print( mul_pool_names[n][t]+"\t" );
    		}
    		System.out.println();
    	}
    	
    	for (int n=0; n< num_run ;n++ ) {
    		int real_num_threads= 0;
    		for (int t=0; t< num_threads ;t++ ) {
    			if (! mul_pool_names[n][t].equals("NULL")) {
    				real_num_threads++;
    			}
    		}
    		System.out.println("Running L0L1 regression for " +Integer.toString(level)+
            		" region "+ Integer.toString(region) +"; "+ 
    		Integer.toString(real_num_threads)+" pools:" );
    		for (int t=0; t<real_num_threads;t ++ ) {
    			try {  
		            String shpath=rscript_path+ " "+ rscript_prefix+"_pool_"+ mul_pool_names[n][t]+".R";  
		            Process ps = Runtime.getRuntime().exec(shpath, null); 
		            System.out.print( mul_pool_names[n][t] +"\t");
    			} catch (Exception e) {  
	                e.printStackTrace();  
    			}
    		}
    		System.out.println( );
    		int iter_time= 1000*60*10/20;
	    	int maximum_iter_time= 2*60*12; 
	    	int time_count =0;
	    	while (time_count< maximum_iter_time) {
	    		time_count++;
	    		int file_count=0;
	    		Thread.sleep(iter_time );
	    		for (int t=0; t<real_num_threads;t ++ ) {
	    			String outfile= outfile_prefix+ mul_pool_names[n][t]+"/regression_level_" + Integer.toString(level)+
					"_region_"+ Integer.toString(region)+ ".regression_out";
	    			File file = new File(outfile);
	    			File file_log = new File(outfile+".log");
	    			if ((file.exists())   || (file_log.exists())) {
	    				file_count++;
	    			}
	    		}
	        	System.out.println("Of "+ Integer.toString(real_num_threads)+" pools above, "+Integer.toString(file_count)+"/"+
	                		Integer.toString( real_num_threads )+
	                		" are finished:\t"+ dtf.format(LocalDateTime.now()));
	    		if (file_count==real_num_threads ) {
	    			break;
	    		}
	    	}
    	}
    		
    		
    		
//	    	for (int p=0; p<pool_names.length;p++ ) {
//		    	try {  
//		            String shpath=rscript_path+ " "+ rscript_prefix+"_pool_"+ pool_names[p]+".R";  
//		            Process ps = Runtime.getRuntime().exec(shpath, null); 
//		            System.out.println(rscript_prefix+"_pool_"+ pool_names[p]+".R");
		                
	//	                ps.waitFor();  
	//	                BufferedReader br = new BufferedReader(new InputStreamReader(ps.getInputStream()));  
	//	                StringBuffer sb = new StringBuffer();  
	//	                String line;  
	//	                while ((line = br.readLine()) != null) {  
	//	                    sb.append(line).append("\n");  
	//	                }  
		//                String result = sb.toString();  
		//                System.out.println(result);               
		//                if (result.contains("FINISH")){
		//                	System.out.println("regional L0L1 regression for level "+ Integer.toString(level)+
		//                		" region "+ Integer.toString(i+1)  +" finished:\t"+ dtf.format(LocalDateTime.now()));
		//                }
	//	                br.close();
//		            }   catch (Exception e) {  
//		                e.printStackTrace();  
//		        }
//	    	}
//	
//	    	
//	    	int iter_time= 1000*60*10/20;
//	    	int maximum_iter_time= 6*24;
//	    	int time_count =0;
//	    	
//	    	while (time_count< maximum_iter_time) {
//	    		int file_count=0;
//	    		Thread.sleep(iter_time );
//	    		for (int p=0; p<pool_names.length;p++ ) {
//	    			String outfile= outfile_prefix+ pool_names[p]+"/regression_level_" + Integer.toString(level)+
//					"_region_"+ Integer.toString(region)+ ".regression_out";
//	    			File file = new File(outfile);
//	    			if (file.exists()) {
//	    				file_count++;
//	    			}
//	    		}
//	        	System.out.println("regional L0L1 regression for level "+ Integer.toString(level)+
//	                		" region "+ Integer.toString(region)  + ":\t"+Integer.toString(file_count)+"/"+
//	                		Integer.toString(pool_names.length)+
//	                		" finished:\t"+ dtf.format(LocalDateTime.now()));
//	    		if (file_count==pool_names.length ) {
//	    			return;
//	    		}
//	    	}
    	
	    	
    }


    /**
     *  TODO: (old) figure out how to initialize weights based on the relative importance of 1,
     *  H and H^H.
     *
     *  @param weights
     */
    public void setup_weights(double[] weights) {
        this.maf_weights = new double[this.num_loci];
        this.ld_weights = new double[this.num_loci][this.num_loci];
        this.sum1_weight = weights[0]; // TODO: (old) think about how to set the weight
        double maf_weight = weights[1]; // TODO: (old) think about how to set these weights
        double ld_weight = weights[2]; // TODO: (old) think about how to set these weights
        for (int loc = 0; loc < this.num_loci; loc++) {
            this.maf_weights[loc] = maf_weight;
        }

        for (int loc1 = 0; loc1 < this.num_loci; loc1++) {
            for (int loc2 = loc1 + 1; loc2 < this.num_loci; loc2++) {
                this.ld_weights[loc1][loc2] = ld_weight;
            }
        }
    }


    /**
     *  Set up a HashMap for search of indexes in the locations array when processing VEF files.
     */
    public void setup_site_locations2index_map() {
        this.site_locations2index = new HashMap<Integer, Integer>();
        for (int loc = 0; loc < this.num_loci; loc++) {
            this.site_locations2index.put(this.potential_haps.locusInfo[loc].start_loc, loc);
        }
    }

    /**
     *  The three functions will calculate ld_matrix differently:
     * 	  (1) calcualte_LD_matrix_single_pool()
     *    (2) calcualte_LD_matrix_multi_pool_reads()
     *    (3) calcualte_LD_matrix_multi_pool_stat()
     *
     *  Then they share this same functions of writing to file in libsvm format (this function) and
     *  lasso regression (conduct_regression()).
     */
    public void write_1_H_HH_to_file(String lasso_in_file) {

        // TODO: [LEFTOVER]
        // // ArrayList<Integer> quanlified_hap_index_array = new ArrayList<Integer>();
        // for (int h = 0; h < this.potential_haps.num_global_hap; h++) {
        //     if (this.potential_haps.global_haps_freq[h] >= this.raw_hap_freq_cutoff) {
        //         quanlified_hap_index_array.add(h);
        //     }
        // }
        // this.potential_haps.num_global_hap = 0;
        // this.quanlified_hap_index = new int[this.potential_haps.num_global_hap];
        // for (int h = 0; h < this.potential_haps.num_global_hap; h++) {
        //     this.quanlified_hap_index[h]=quanlified_hap_index_array.get(h);
        // }

        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(lasso_in_file));

            // Write 1-vector. Format is w_allele\sh1:w\sh2:w...hn:w\n.
            bw.write(Integer.toString(0)+"| "+this.sum1_weight + "");

            // Originally, num_global_hap.
            for (int h = 0; h < this.potential_haps.num_global_hap; h++) {
//                bw.write(" " + (h + 1) + ":" + this.sum1_weight);
                bw.write(" "  + this.sum1_weight);
            }

            bw.write("\n");

            // Write the block of MAF and H.
            for (int loc=0;loc<this.num_loci;loc++) {
                bw.write( Integer.toString(this.index_2_pos_dict.get(loc))+"|"+
                		" "+this.maf_weights[loc] * this.maf[loc] + "");
                for(int h = 0; h < this.potential_haps.num_global_hap; h++) {

                    // TODO: [LEFTOVER]
                    // int hap_index = this.quanlified_hap_index[h];

                    bw.write(" "
//                        + (h + 1)
//                        + ":"
                        + this.maf_weights[loc]

                        // MAF_weights are the relative importance factor assigned to the variant
                        // composition of the potential haplotypes.
                        * this.potential_haps.global_haps[h][loc]);

                }
                bw.write("\n");
            }

            // Write the block of LD and H^H.
            for (int loc1 = 0; loc1 < this.num_loci; loc1++) {
                for (int loc2 = loc1 + 1; loc2 < this.num_loci; loc2++) {
                    // In cases where there is zero sequencing coverage of two sites, there is no
                    // information i.e.: NaN in the 'LD' matrix.
                    if (!Double.isNaN(this.ld_matrix[loc1][loc2]))   {
                    	boolean flag= true;
                    	
                    	
//                    	if (!this.sequencing_technology.equals("10X_linked_reads")) {
//                    		flag =true;
//                    	}
                    	if (flag) {
	                    	double coverage_coefficient= this.lasso_coverage_weight* 
	                    			(double)this.total_counts[loc1][loc2]* 
	                    			this.total_ld_covered_paired_alleles/ this.total_ld_covered_reads_counts;
	                    	double dist_coefficient= 1+ ((double) (loc2-loc1)/ (double)this.num_loci)*
	                    			(this.lasso_distance_max_weight-1);
	                    	double cross_region_weight= 1;
//	                    	if (!this.sequencing_technology.equals("10X_linked_reads")) {
//	                    		if (!this.pos_come_from_region.get(this.index_2_pos_dict.get(loc1)).
//	                    			equals(this.pos_come_from_region.get(this.index_2_pos_dict.get(loc2)))) {
//	                    			cross_region_weight= cross_region_weight* this.lasso_cross_region_weight;
//	                    		}
//	                    	}
	                    	
	//                    	coverage_coefficient=1;
	//                    	dist_coefficient=1;
	                        bw.write(Integer.toString(this.index_2_pos_dict.get(loc1))+":"+ 
	                        		Integer.toString(this.index_2_pos_dict.get(loc2))+"|"
	                        		+" "+this.ld_weights[loc1][loc2] * this.ld_matrix[loc1][loc2]* 
	                        		coverage_coefficient*dist_coefficient *cross_region_weight + "");
	                        for (int h = 0; h < this.potential_haps.num_global_hap; h++) {
	
	                            // TODO: [LEFTOVER]
	                            // int hap_index=this.quanlified_hap_index[h];
	
	                            // TODO: DONE?(old) Lauren: please double-check whether Quan's line
	                            // below is correct!!!!
	                            // Below is the correlation, or co-occurrence of alternate alleles at
	                            // loci 1 and 2 on the same haplotype.
	                        	
	                            bw.write(" "
	//                                + (h + 1)
	//                                + ":"
	                                + this.ld_weights[loc1][loc2]
	                                *cross_region_weight
	                                * dist_coefficient
	                                * coverage_coefficient
	                                * this.potential_haps.global_haps[h][loc1]
	                                * this.potential_haps.global_haps[h][loc2]);
	
	                        }
	                        bw.write("\n");
                    	}
                    }
                }
            }

            bw.close();

        } catch(Exception e) {
            e.printStackTrace();
        }
    }


    /**
     *  copy the MAF from this.potential_haps.inpool_site_freqs
     */
    
    public void get_MAF_single_pool() {
        if (pool_index == -1) {
            System.out.println("ERROR: pool_index==-1; should not calculate single_pool MAF.");
        }

        this.maf = new double[this.num_loci];
        for (int loc = 0; loc < this.num_loci; loc++) {
            this.maf[loc] = this.potential_haps.inpool_site_freqs[loc][this.pool_index];
        }
    }
    
    /**
     *  Calculate global MAF, which is the average of all in-pool MAF at the same site.
     */
    public void calculate_MAF_multi_pool() {
        if (pool_index != -1) {
            System.out.println("ERROR: pool_index!=-1; should not calculate global MAF.");
        }

        this.maf = new double[this.num_loci];
        for (int loc = 0; loc < this.num_loci; loc++) {
            this.maf[loc] = 0;
            for (int pool = 0; pool < this.potential_haps.num_pools; pool++) {
                this.maf[loc] += this.potential_haps.inpool_site_freqs[loc][pool];
            }

            this.maf[loc] = this.maf[loc] / this.potential_haps.num_pools;
        }
    }


    /**
     *  Given any two bi-allelic loci i and j, calculate P(i==1 && j==1).
     *
     *  For single pool, just check the number of reads (including paired-end ones) that cover both
     *  loci. The above probability is then #(1,1)/Sum(#(1,1)+#(0,1)+#(1,0)+#(0,0)).
     *
     *  The string the parse in the VEF files are in the format location=allele; e.g.:
     *
     *  272=0;288=0;583=0;601=1;631=0;  (they are in the 2nd column of the VEF file)
     *
     *  TODO: (old) assign different ld_weights based on different sequencing coverage.
     *
     *  @param vef_file
     *  @return
     */
    public double[][] calcualte_LD_matrix_single_pool(String vef_file) {

        // TODO: [LEFTOVER]
        // if (pool_index == -1) {
        //     System.out.println("ERROR: pool_index==-1; should not calculate single_pool LD.");
        // }
        //
        // // The array this.ld_matrix will initially store counts of (1,1); the array total_counts
        // // store all counts. Finally this.ld_matrix/total_counts will be the outcome.
        // // this.ld_matrix[loc1][loc2] will be NaN if total_counts[loc1][loc2]==0;
        // this.setup_site_locations2index_map();

        // this.num_loci is already synced up with the input HapConfig object.
    	
    	
        double[][] local_matrix = new double[this.num_loci][this.num_loci];
        HashSet<String> pos_pairs_set = new HashSet<String>();
        this.total_ld_covered_paired_alleles=0;
        this.total_counts = new int[this.num_loci][this.num_loci];
        try {
            BufferedReader br = new BufferedReader(new FileReader(vef_file));
            String line = br.readLine();
            while (line!=null) {
                String[] location_alleles = line.split("\t")[1].split(";");
                if (location_alleles.length < 2) { // no multiple sites in this read
                    line = br.readLine();
                    continue;
                }

                // The search may fail, so we don't know the length of the array.
                ArrayList<Integer> indexes = new ArrayList<Integer>();
                ArrayList<String> alleles = new ArrayList<String>();
                ArrayList<Integer> positions = new ArrayList<Integer>();
                for (int i = 0; i < location_alleles.length; i++) {
                    String[] tmp_loc_allele = location_alleles[i].split("=");
                    int the_location = Integer.parseInt(tmp_loc_allele[0]);

                    // TODO: [LEFTOVER]
                    // setup_site_locations2index_map<POS_ON_REF,POS_IND_HAPCONFIG>

                    // If we do find the site in question...
                    if (this.site_locations2index.containsKey(the_location)) {
                        indexes.add(this.site_locations2index.get(the_location));
                        alleles.add(tmp_loc_allele[1]);
                        positions.add(the_location);
                    }
                }
                
                

                int num_avail_locs=indexes.size();
                for (int site1 = 0; site1 < num_avail_locs; site1++) {
                    for (int site2 = site1 + 1; site2 < num_avail_locs; site2++) {
                        int index1 = indexes.get(site1);
                        int index2 = indexes.get(site2);
//                        System.out.println(positions.get(site1) +"\t"+
//                        		positions.get(site2));
//                        System.out.println(this.pos_come_from_region.get(positions.get(site1)) +"\t"+
//                        		this.pos_come_from_region.get(positions.get(site2)));
                       if ((!this.sequencing_technology.equals("10X_linked_reads"))  ||
                    		   ( !this.pos_come_from_region.get(positions.get(site1)).
                    		   equals( this.pos_come_from_region.get(positions.get(site2)) ))  || 
                       (this.num_regions_from==1)  || 
                       (this.l0l1_regression_level==1)) {
//                    	   System.out.println(this.pos_come_from_region.get(positions.get(site1)) +"\t"+
//                           		this.pos_come_from_region.get(positions.get(site2)));
	                        if (index1 > index2) { // ld_matrix only record loc1 < loc2 entries
	                            int tmp = index1;
	                            index1 = index2;
	                            index2 = tmp;
	                        }
	                        //Chen:
	                        this.total_counts[index1][index2]++;
	                        this.total_ld_covered_reads_counts =
	                        		this.total_ld_covered_reads_counts+1;
	                        String tmp = Integer.toString(index1)+":"+ Integer.toString(index2);
	                        if (!pos_pairs_set.contains(tmp)) {
	                        	pos_pairs_set.add(tmp);
	                        	this.total_ld_covered_paired_alleles= 
	                        			this.total_ld_covered_paired_alleles+1;
	                        }
	                        // Both are "1".
	                        if (alleles.get(site1).equals("1") && alleles.get(site2).equals("1")) {
	                            local_matrix[index1][index2]++;  
	                        }
                       }
                    }
                }

                line=br.readLine();
            }

            br.close();

            // Divide by the total number of reads covering both sites.
            for (int loc1 = 0; loc1 < this.num_loci; loc1++) {
                for (int loc2 = 0; loc2 < this.num_loci; loc2++) {
                    if (this.total_counts[loc1][loc2] != 0) {
                        local_matrix[loc1][loc2] = local_matrix[loc1][loc2]
                            / this.total_counts[loc1][loc2];
                        
                    } else {
                        local_matrix[loc1][loc2] = Double.NaN;
                    }
                }
            }

        } catch(Exception e) {
            e.printStackTrace();
        }

        return local_matrix;
    }


    /**
     *  Given any two bi-allelic loci i and j, calculate P(i==1 && j==1).
     *
     *  For many pools, aggregate the LD_matrix from individual pools (that were estimated by
     *  reads-counting).
     *
     *  Before entering the function, one has to ENSURE that the order of pools in "potential_haps"
     *  is the same as the vef_files[]!!!
     *
     *  @param vef_file
     */
    public void calcualte_LD_matrix_multi_pool_reads(String[] vef_file) {

        // TODO: [LEFTOVER]
        // Regularization[] single_pools=new Regularization[this.potential_haps.num_pools];

        this.setup_site_locations2index_map();
        int[][] pool_counter = new int[this.num_loci][this.num_loci];
        this.ld_matrix = new double[this.num_loci][this.num_loci];
        for (int pool = 0; pool < this.potential_haps.num_pools; pool++) {

            // TODO: [LEFTOVER]
            // single_pools[pool] = new Regularization(
            //     pool_index,
            //     Double.NaN,
            //     potential_haps,
            //     raw_hap_freq_cutoff,
            //     null,
            //     null);
            //
            // // This is unnecessary; but put here is case new versions of this function that
            // // involves MAF may be developed.
            // single_pools[pool].get_MAF_single_pool();

            double[][] local_ld_matrix = calcualte_LD_matrix_single_pool(vef_file[pool]);

            // TODO: [LEFTOVER]
            // // And then add it to the global ld_matrix.
            // for (int loc1 = 0; loc1 < this.num_loci; loc1++) {
            //     for (int loc2 = loc1 + 1; loc2 < this.num_loci; loc2++) {
            //         System.out.println(this.potential_haps.locusInfo[loc1].start_loc
            //             + "\t"
            //             + this.potential_haps.locusInfo[loc2].start_loc
            //             + "\t"
            //             + local_ld_matrix[loc1][loc2]);
            //     }
            //
            // // Calculate the ld_matrix in the single pool;
            // single_pools[pool].calcualte_LD_matrix_single_pool(vef_file[pool]);

            // And then add it to the global ld_matrix.
            for (int loc1 = 0; loc1 < this.num_loci; loc1++) {
                for (int loc2 = loc1 + 1; loc2 < this.num_loci; loc2++) {
                    if (local_ld_matrix[loc1][loc2] != Double.NaN) {
                        this.ld_matrix[loc1][loc2] += local_ld_matrix[loc1][loc2];

                        // TODO: [LEFTOVER]
                        // single_pools[pool].ld_matrix

                        pool_counter[loc1][loc2]++;
                    }
                }
            }
        }
        for (int loc1 = 0; loc1 < this.num_loci; loc1++) { // normalize the global ld_matrix
            for (int loc2 = loc1 + 1; loc2 < this.num_loci; loc2++) {
                this.ld_matrix[loc1][loc2] = this.ld_matrix[loc1][loc2]
                    / (double) pool_counter[loc1][loc2];

                // TODO: [LEFTOVER]
                // this.potential_haps.num_pools;

            }
        }
    }

    /**
     *  Given any two bi-allelic loci i and j, calculate P(i==1 && j==1).
     *
     *  For many pools, calculate the co-variance between each loci using the formula in Zhang H,
     *  Yang HC, and Yang YN Bioinformatics 2008 paper.
     *
     *  variance–covariance matrix \Sigma = (1/#pool)(SUM(A_i x A_i' - \omiga x \omiga'), which
     *  essentially is COV(X,Y)=E(XY)-E(x)E(y).
     */
    public void calcualte_LD_matrix_multi_pool_stat() {
        if (pool_index != -1) {
            System.out.println("ERROR: pool_index!=-1; should not calculate global LD.");
        }

        this.ld_matrix = new double[this.num_loci][this.num_loci];
        double[][] inpool_site_freqs = this.potential_haps.inpool_site_freqs;
        for (int loc1 = 0; loc1 < this.num_loci; loc1++) {
            for (int loc2 = loc1 + 1; loc2 < this.num_loci; loc2++) {
                // TODO: Lauren: check whether the following four lines that provide the
                // calculation of \Sigma_0 (in Zhang et al 2008 Bioinfo) is correct for Quan!!!
                for (int pool = 0; pool < this.potential_haps.num_pools; pool++) {
                    // SUM(XY)
                    this.ld_matrix[loc1][loc2] += inpool_site_freqs[loc1][pool]
                        * inpool_site_freqs[loc1][pool];
                }

                // E(XY) = SUM(XY) / n
                this.ld_matrix[loc1][loc2] = this.ld_matrix[loc1][loc2]
                    / this.potential_haps.num_pools;

                // COV(XY) = E(XY) - E(X)E(Y)
                this.ld_matrix[loc1][loc2] = this.ld_matrix[loc1][loc2]
                    - this.maf[loc1] * this.maf[loc2];
            }
        }
    }


    /**
     *
     * @throws FileNotFoundException
     */
    public void conduct_regression(String lasso_in, String lasso_out) throws FileNotFoundException {
    	/*     turn off spark
        // Initiate the Spark session.
        PrintStream originalOut = System.out;
        PrintStream originalErr = System.err;

        File file = new File(lasso_out);
        PrintStream fileOut = new PrintStream(file);
        PrintStream fileErr = new PrintStream(file);
        System.setOut(fileOut);
        System.setErr(fileErr);

        SparkSession spark = SparkSession.builder()
            .appName("sparkAnalysis")
            .master("local[*]")
            .config("spark.driver.memory", this.memory_usage)
            .getOrCreate();

        // Load training data.
        Dataset<Row> training = spark.read().format("libsvm").load(lasso_in);

        // Create the regression object with L1 regularization and no intercept.
  
        LinearRegression lr = new LinearRegression()
            .setMaxIter(10)
            .setRegParam(this.lambda)
            .setElasticNetParam(1) 
            .setFitIntercept(false);

        // Fit the model.
        LinearRegressionModel lrModel = lr.fit(training);

        // TODO: [LEFTOVER]
        // System.out.println("LinearRegression DONE");

        this.out_hap_freqs = lrModel.coefficients().toArray();

        // Below are some properties of the LASSO regression for debug purposes.
        // Summarize the model over the training set and print out some metrics.
        fileOut.close();
        fileErr.close();
        System.setOut(originalOut);
        System.setErr(originalErr);
        LinearRegressionTrainingSummary trainingSummary = lrModel.summary();
        System.out.println("numIterations: " + trainingSummary.totalIterations());
        System.out.println("objectiveHistory: "
            + Vectors.dense(trainingSummary.objectiveHistory()));

        trainingSummary.residuals().show();
        System.out.println("RMSE: " + trainingSummary.rootMeanSquaredError());
        System.out.println("r2: " + trainingSummary.r2());
        //this.r2 = trainingSummary.r2();
         * 
         */
    }


    /**
     *  The main pipeline the links other functions.
     *
     *  There are three branches:
     *  (1) Single pool. (Has to be read-based ld_matrix calculation.)
     *  (2) Multiple pools, statistics-based ld_matrix calculation.
     *  (3) Mutiple pools, read-based ld_matrix calculation.
     *
     *  @param vef_file
     *  @param vef_files
     *  @param weights
     *  @throws FileNotFoundException
     */
    
    public void estimate_frequencies_lasso(
        String vef_file,
        String[] vef_files,
        double[] weights) throws FileNotFoundException {
    	
        // First calculate Ys (MAF and LD) from pooled data.
        if (vef_file != null && vef_files == null) { // there is a vef_file, indicating single pool
            if (this.pool_index == -1) {
                System.out.println("ERROR: pool_index==-1 but we do have a VEF file");
            }
            this.get_MAF_single_pool();
            this.setup_site_locations2index_map();
            this.ld_matrix = calcualte_LD_matrix_single_pool(vef_file);
        } else if (vef_file == null) { // multiple pools.
            if (this.pool_index != -1) {
                System.out.println("ERROR: pool_index!=-1 but we don't have a VEF file");
            }
            this.calculate_MAF_multi_pool();
            if(vef_files!=null) { // multiple VEF files supplied, indicating read-based
                this.calcualte_LD_matrix_multi_pool_reads(vef_files);

            } else {  // no VEF files, indicating statistic-based
                this.calcualte_LD_matrix_multi_pool_stat();
            }

            // TODO: (old) in the future, we may have a version use both stat- and read-based
            // ld_matrix and a parameter to structure their relative weights.

        } else {
            System.out.println("ERROR: Combination of parameters wrong!");
        }

        // Set up weights after calculating Ys.
        // TODO: (old) note that setup_weights() should be implemented differently based on the
        // type of calculating ld_matrix in the future.
        this.setup_weights(weights);

        // Prepare file in libsvm format ready for Spark LASSO.
        this.write_1_H_HH_to_file(this.prefix + ".regression_in");

        // Run LASSO to get the output hap frequencies.
//        this.conduct_regression(this.prefix + ".lasso_in", this.prefix + ".lasso_out");
    }
    
    
    
    public void estimate_frequencies_lasso(
            String vef_file,
            String[] vef_files,
            double[] weights,
            String sequencing_tech, 
            HashMap<Integer,String>  pos_come_from,
            int  num_regions_come_from, 
            HashMap<Integer,Integer>  index_2_pos, 
            int l0l1_level)
            		throws FileNotFoundException {
    		this.sequencing_technology= sequencing_tech;
    		this.pos_come_from_region = new HashMap<Integer,String> ();
    		this.index_2_pos_dict = new HashMap<Integer, Integer> ();
    		this.pos_come_from_region= pos_come_from;
    		this.num_regions_from= num_regions_come_from;
    		this.index_2_pos_dict= index_2_pos;
    		this.l0l1_regression_level= l0l1_level+1;
        	if (sequencing_tech.equals("10X_linked_reads")) {
        		this.lasso_distance_max_weight=1;
        	}
        	
            // First calculate Ys (MAF and LD) from pooled data.
            if (vef_file != null && vef_files == null) { // there is a vef_file, indicating single pool
                if (this.pool_index == -1) {
                    System.out.println("ERROR: pool_index==-1 but we do have a VEF file");
                }
                this.get_MAF_single_pool();
                this.setup_site_locations2index_map();
                this.ld_matrix = calcualte_LD_matrix_single_pool(vef_file);
            } else if (vef_file == null) { // multiple pools.
                if (this.pool_index != -1) {
                    System.out.println("ERROR: pool_index!=-1 but we don't have a VEF file");
                }
                this.calculate_MAF_multi_pool();
                if(vef_files!=null) { // multiple VEF files supplied, indicating read-based
                    this.calcualte_LD_matrix_multi_pool_reads(vef_files);

                } else {  // no VEF files, indicating statistic-based
                    this.calcualte_LD_matrix_multi_pool_stat();
                }

                // TODO: (old) in the future, we may have a version use both stat- and read-based
                // ld_matrix and a parameter to structure their relative weights.

            } else {
                System.out.println("ERROR: Combination of parameters wrong!");
            }

            // Set up weights after calculating Ys.
            // TODO: (old) note that setup_weights() should be implemented differently based on the
            // type of calculating ld_matrix in the future.
            this.setup_weights(weights);

            // Prepare file in libsvm format ready for Spark LASSO.
            this.write_1_H_HH_to_file(this.prefix + ".regression_in");

            // Run LASSO to get the output hap frequencies.
//            this.conduct_regression(this.prefix + ".lasso_in", this.prefix + ".lasso_out");
        }

    public HapConfig hapOut(String[] pool_IDs) {
    	//String[] single_pool_IDs={pool_name};
        HapConfig raw_final_haps = new HapConfig(
            this.potential_haps.global_haps_string,
            this.out_hap_freqs,
            null,
            this.potential_haps.inpool_site_freqs,
            this.potential_haps.locusInfo,
            1,
            this.potential_haps.hap_IDs,
            pool_IDs,
            this.potential_haps.est_ind_pool);
        // This gets rid of all haplotypes that have 0-frequency.
        
        return raw_final_haps.filter(this.raw_hap_freq_cutoff);

        // TODO: [LEFTOVER]
        // final_haps.write_global_file_string(vef_file + "_2.inter_freq_vars.txt", false);

    }
}
