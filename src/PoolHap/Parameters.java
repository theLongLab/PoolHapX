package PoolHap;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashSet;
import java.util.Properties;

public class Parameters {
    // TODO: [ReconEP]:: add any new parameters as needed.

    public static String[] supported_functions_array = {"format", "gc", "aem", "lasso"};
    public String function; // which module to run: format, gc, aem or lasso.

    /*
     *  General parameter set.
     */    

    public String project_name;
        // input files: SAM files and VCF files. Added by Quan Long 2019-07-01
    public String input_dir;
    public String inter_dir; // intermediate directory, including the following files:
    public String out_dir; // output directory
    public String gold_dir;
    public int num_pos_window;
    public int num_gap_window; 
    public int num_pos_job;
    //  Divide and conquer parameter set.
    public double gap_inpool_cutoff; // a ratio
    public double gap_all_pool_cutoff; // a ratio
//    public double gap_support_step; // a ratio
    public int min_level_I_region_size;
    public int max_level_I_region_size;
//    public int min_level_I_last_size; // Quan removed it July 27, 2019.
    public int min_level_II_region_size;
    public int max_level_II_region_size;
    public int est_ind_pool;
    
    public int virtual_cov_link_gc;  // how many total "coverage" in the pool. 
                                        //The larger, the higher resolution for rare haplotypes of 
                                        //frequency 1.0/virtual_cov_link_gc.
    
    

    // Approximate expectation-maximization parameter set.
    public int aem_max_iteration;
    public double aem_epsilon;
    public double aem_zero_cutoff;
//    public double final_cutoff;
//    public int adhoc_freq_cutoff;
    public double aem_regional_cross_pool_freq_cutoff;
    public int aem_hapset_size_max;
    public int aem_hapset_size_min;
    public int aem_maximum_level;
    public double if_0_0;
    public double if_Denominator_0;
   
    // LASSO parameters    
    // global
    public double lasso_full_hap_freq_cutoff;
   
    public double lasso_global_lambda;
    public double[] lasso_weights;
    public String lasso_global_memory;
    
    // regional
    public double lasso_regional_lambda;
    public String lasso_regional_memory;
    public double lasso_regional_cutoff;
    public double lasso_regional_cross_pool_cutoff;
    public int lasso_hapset_size_max;
    public int lasso_hapset_size_min;
    	
    public int bfs_mismatch_tolerance;
    public int level_III_IV_region_mismatch_tolerance;
    public int level_V_VI_region_mismatch_tolerance;
    public int level_VII_VIII_region_mismatch_tolerance;
    public double lasso_coverage_weight;
    public double lasso_distance_max_weight;
    public double hc_similarity_cutoff;
    public double mcc_freq_cutoff;
    
    public String rscript_path;
    public double regression_gamma_min;
    public double regression_gamma_max;
    public int regression_hapset_size_max;
    public int ngammas;
    public int regression_link_mismatch_tolerance;
    public int regression_hapset_size_min;
    public int regression_maximum_regions;
    public int regression_maximum_selected_haplotypes;
    public String sequencing_technology; 
    public int num_threads;
    public String species;
    
    
    //public double hapset_size_rand;
    //public int min_num_hap_regiobal_lasso;
    
    public Parameters(String propFilePath) throws IOException {
            
//        HashSet<String> supported_functions = new HashSet<String>();
//        for (int k = 0; k < supported_functions_array.length; k++) {
//             supported_functions.add(supported_functions_array[k]);
//        }

        Properties prop = new Properties(); // properties object from properties file
        InputStream is = new FileInputStream(propFilePath);
        prop.load(is);
        
//        this.function = prop.getProperty("Function");
//        if (!supported_functions.contains(this.function)) {
//            System.out.println("Function "+this.function+" is not supported. A typo?");
//            System.exit(0);
//        }
        
        this.project_name= prop.getProperty("Proj_Name");
        this.input_dir= prop.getProperty("Input_Dir")+"/";
        this.inter_dir = prop.getProperty("Intermediate_Dir")+"/";
        this.out_dir = prop.getProperty("Output_Dir")+"/";
        this.gold_dir= prop.getProperty("Gold_Dir")+"/";
        this.num_pos_window = Integer.parseInt(prop.getProperty("Num_Pos_Window"));
        this.num_gap_window = Integer.parseInt(prop.getProperty("Num_Gap_Window"));
        this.num_pos_job= 1000; 
        // Divide and conquer 
        
        this.gap_inpool_cutoff = Double.parseDouble(
            prop.getProperty("In-pool_Gap_Support_Min"));
        this.gap_all_pool_cutoff =  Double.parseDouble(
            prop.getProperty("All-pool_Gap_Support_Min"));
        this.min_level_I_region_size = Integer.parseInt(
            prop.getProperty("Level_1_Region_Size_Min"));
        this.max_level_I_region_size = Integer.parseInt(
            prop.getProperty("Level_1_Region_Size_Max"));
 //       this.min_level_I_last_size = Integer.parseInt(
 //           prop.getProperty("Level_1_Last_Region_Min"));
        this.min_level_II_region_size = Integer.parseInt(
            prop.getProperty("Level_1T_Region_Size_Min"));
        this.max_level_II_region_size = Integer.parseInt(
            prop.getProperty("Level_1T_Region_Size_Max"));
        this.est_ind_pool = Integer.parseInt(prop.getProperty("Est_Ind_PerPool"));
        this.if_0_0 = Double.parseDouble
        		(prop.getProperty("IF_0_0")); 
        
        this.if_Denominator_0 = Double.parseDouble
        		(prop.getProperty("IF_Denominator_0")); 
        
        this.level_III_IV_region_mismatch_tolerance = Integer.parseInt(
                prop.getProperty("Level_1_2_Region_Mismatch_Tolerance"));
        this.level_V_VI_region_mismatch_tolerance = Integer.parseInt(
                prop.getProperty("Level_2_3_Region_Mismatch_Tolerance"));
        this.level_VII_VIII_region_mismatch_tolerance = Integer.parseInt(
                prop.getProperty("Level_3_4_Region_Mismatch_Tolerance"));
        this.bfs_mismatch_tolerance = Integer.parseInt(
                prop.getProperty("BFS_Mismatch_Tolerance"));
        
        // GC Link regions:
//        this.virtual_cov_link_gc = Integer.parseInt(prop.getProperty("Virtual_Coverage_Link_GC"));
        this.hc_similarity_cutoff = 0.95;
        
        // LASSO
//        this.lasso_global_lambda = Double.parseDouble(prop.getProperty("LASSO_Global_Lambda_Penalty"));
//        this.lasso_regional_lambda = Double.parseDouble(prop.getProperty("LASSO_Reginal_Lambda_Penalty"));
//        this.lasso_global_memory = prop.getProperty("LASSO_Global_Memoery");
//        this.lasso_regional_memory = prop.getProperty("LASSO_Regional_Memoery");
        this.lasso_weights = new double[] {
            Double.parseDouble(prop.getProperty("Regression_One_Vector_Weight")),
            Double.parseDouble(prop.getProperty("Regression_Hap_MAF_Weight")),
            Double.parseDouble(prop.getProperty("Regression_Hap_LD_Weight"))};
        
        this.lasso_coverage_weight = Double.parseDouble(
                prop.getProperty("Regression_Coverage_Weight"));
        this.lasso_distance_max_weight = Double.parseDouble(
                prop.getProperty("Regression_Distance_Max_Weight"));

//        this.min_r2 = Double.parseDouble(prop.getProperty("Minimum_R2_Fit"));
//        this.lasso_penalty_step = Double.parseDouble(prop.getProperty("Penalty_Step_Size"));
//        this.lasso_regional_cross_pool_cutoff = Double.parseDouble(
//            prop.getProperty("Regional_Global_Freq_Min"));
//        this.lasso_full_hap_freq_cutoff = Double.parseDouble(
//            prop.getProperty("LASSO_Full_Length_Inpool_Freq_Min"));
//        this.lasso_regional_cross_pool_cutoff=Double.parseDouble(
//            prop.getProperty("LASSO_Regional_Cross_Pool_Freq_Min"));
//        this.lasso_hapset_size_max = Integer.parseInt(
//            prop.getProperty("LASSO_Regional_HapSetSize_Max"));
//
//        this.lasso_hapset_size_min = Integer.parseInt(
//            prop.getProperty("LASSO_Regional_HapSetSize_Min"));

//        this.hapset_size_rand = Double.parseDouble(prop.getProperty("DC_HapSetSize_Rand"));
        
        /*
         *  Extract parameters to AEM parameter object variables from properties object.
         */
        this.aem_max_iteration = Integer.parseInt(prop.getProperty("AEM_Iterations_Max"));
        this.est_ind_pool = Integer.parseInt(prop.getProperty("Est_Ind_PerPool"));
        this.aem_epsilon = Double.parseDouble(prop.getProperty("AEM_Convergence_Cutoff"));
        this.aem_zero_cutoff = Double.parseDouble(prop.getProperty("AEM_Zero_Cutoff"));
        this.aem_maximum_level = Integer.parseInt(prop.getProperty("AEM_Maximum_Level"))*2-1;
        this.aem_regional_cross_pool_freq_cutoff = Double.parseDouble(
            //prop.getProperty("Regional_Global_Freq_Min"));
            prop.getProperty("AEM_Regional_Cross_Pool_Freq_Cutoff"));
        this.aem_hapset_size_max = Integer.parseInt(
            prop.getProperty("AEM_Regional_HapSetSize_Max"));
        this.aem_hapset_size_min = Integer.parseInt(
            prop.getProperty("AEM_Regional_HapSetSize_Min"));
        this.mcc_freq_cutoff = 0.01;
        
        this.regression_gamma_max=Double.parseDouble
		(prop.getProperty("Regression_Gamma_Max"));
        this.regression_gamma_min=Double.parseDouble
		(prop.getProperty("Regression_Gamma_Min"));
        this.regression_hapset_size_max= 20;
        this.rscript_path= prop.getProperty("Rscript_path");
        this.ngammas= Integer.parseInt(
                prop.getProperty("Regression_n_Gamma"));
        this.regression_link_mismatch_tolerance= Integer.parseInt(
                prop.getProperty("Regression_Mismatch_Tolerance"));
        this.regression_hapset_size_min=  this.aem_hapset_size_max/2;
//        Integer.parseInt(prop.getProperty("Regression_Regional_HapSetSize_Min"));
       this.regression_maximum_regions=
                Integer.parseInt(prop.getProperty("Regression_Maximum_Regions"));
       
       this.regression_maximum_selected_haplotypes= 20;
       this.sequencing_technology= prop.getProperty("Sequencing_Technology");
       this.num_threads= Integer.parseInt(prop.getProperty("Number_Threads"));
//       this.species= prop.getProperty("Species");
//       
//       if (this.species.equals("bacteria")) {
//    	   this.lasso_weights[0] = this.lasso_weights[0] *2.5;
//       }
       this.species= "virus";
 //      this.aem_hapset_size_rand = Integer.parseInt(prop.getProperty("AEM_HapSetSize_Rand"));
 //      this.adhoc_freq_cutoff = Integer.parseInt(prop.getProperty("Adhoc_Freq_Cutoff"));
        is.close(); // close input stream         
              
        
     }

}


    
       

       
                
  
