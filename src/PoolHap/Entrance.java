package PoolHap;

import java.util.ArrayList;
import java.io.InputStreamReader;
import MiscFunctions.BAMFormatterGATK;
import java.time.temporal.TemporalAccessor;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.HashSet;
import java.util.Iterator;
import java.io.Reader;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.Writer;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Arrays;
import java.io.File;
import java.util.HashMap;

public class Entrance
{
    public static String[] suffixes;
    public static int num_pools;
    public static HashMap<String, Integer> name2index;
    public static String[] names_array;
    
    static {
        Entrance.suffixes = new String[] { "sam", "vcf", "vef", "gcf" };
        Entrance.num_pools = -1;
        Entrance.name2index = new HashMap<String, Integer>();
    }
    
    public static String[] get_filepaths(final String name_file, final String folder, final String suffix, final boolean write_name_file) {
        String[] filepaths = null;
        try {
            final File[] files = new File(folder).listFiles();
            if (write_name_file) {
                filepaths = new String[files.length];
                for (int k = 0; k < files.length; ++k) {
                    filepaths[k] = files[k].toString();
                }
                Arrays.sort(filepaths);
                Entrance.num_pools = files.length;
                Entrance.names_array = new String[Entrance.num_pools];
                for (int index = 0; index < Entrance.num_pools; ++index) {
                    final String[] separat_dirs = filepaths[index].split("/");
                    final String[] separated = separat_dirs[separat_dirs.length - 1].split("\\.");
                    if (!separated[separated.length - 1].equals(suffix)) {
                        System.out.println("ERROR: The suffix of " + filepaths[index] + " is not " + suffix);
                        System.exit(0);
                    }
                    else {
                        String name = separated[0];
                        if (separated.length >= 3) {
                            for (int i = 1; i < separated.length - 1; ++i) {
                                name = String.valueOf(name) + "." + separated[i];
                            }
                        }
                        Entrance.name2index.put(name, index);
                        Entrance.names_array[index] = name;
                    }
                }
                final BufferedWriter bw = new BufferedWriter(new FileWriter(name_file));
                for (int p = 0; p < Entrance.num_pools; ++p) {
                    bw.write(String.valueOf(Entrance.names_array[p]) + "\n");
                }
                bw.close();
            }
            else {
                final BufferedReader br = new BufferedReader(new FileReader(name_file));
                String line = br.readLine();
                int p_index = 0;
                while (line != null) {
                    Entrance.name2index.put(line, p_index++);
                    line = br.readLine();
                }
                br.close();
                Entrance.num_pools = p_index;
                Entrance.names_array = new String[Entrance.num_pools];
                for (final String name : Entrance.name2index.keySet()) {
                    Entrance.names_array[Entrance.name2index.get(name)] = name;
                }
                if (Entrance.num_pools != files.length) {
                    System.out.println("ERROR: No. of files in " + folder + " is " + files.length + ", NOT consistent to the No. of lines in the " + name_file + ": " + Entrance.num_pools);
                    System.exit(0);
                }
                filepaths = new String[files.length];
                final String[] observed_files = new String[files.length];
                for (int j = 0; j < files.length; ++j) {
                    observed_files[j] = files[j].toString();
                }
                for (int index2 = 0; index2 < Entrance.num_pools; ++index2) {
                    final String[] separated2 = observed_files[index2].split("/")[observed_files[index2].split("/").length - 1].split("\\.");
                    if (!separated2[separated2.length - 1].equals(suffix)) {
                        System.out.println("ERROR: The suffix of " + observed_files[index2] + " is not " + suffix);
                        System.exit(0);
                    }
                    else {
                        String name2 = separated2[0];
                        if (separated2.length >= 3) {
                            for (int l = 1; l < separated2.length - 1; ++l) {
                                name2 = String.valueOf(name2) + "." + separated2[l];
                            }
                        }
                        if (!Entrance.name2index.containsKey(name2)) {
                            System.out.println("ERROR: The name of " + observed_files[index2] + " is not in " + "existing name2index table specified by previsous file sets");
                            System.exit(0);
                        }
                        else {
                            final int prespecified_index = Entrance.name2index.get(name2);
                            filepaths[prespecified_index] = observed_files[index2];
                        }
                    }
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return filepaths;
    }
    
    public static void main(final String[] args) throws Exception {
    	
    	
    	
        
        
        final String[] supported_functions_array = { "non_perfect", "format", "gc", "aem", "lasso", "split", "clustering", "evaluate", "analysis", "l0l1", "Comparison_bacteria", "Comparison_metagenomics", "Comparison_human", "hippo", "complete_analysis", "script" };
        final HashSet<String> supported_functions = new HashSet<String>();
        for (int k = 0; k < supported_functions_array.length; ++k) {
            supported_functions.add(supported_functions_array[k]);
        }
        final String function = args[0];
        if (!supported_functions.contains(function)) {
            System.out.println("Function " + function + " is not supported. A typo?");
            System.exit(0); 
        }
        if (function.equals("script")) {
            final String config_path = args[1];
            final Script sc = new Script(config_path);
            sc.ComparisonSim();
            System.out.println("Now, you can run the file in the cmd folder!");
            System.exit(0);
        }
        final String parameter_file = args[1];
        final Parameters gp = new Parameters(parameter_file);
        final String gs_var_pos = String.valueOf(gp.inter_dir) + gp.project_name + "_vars.intra_freq.txt";
        final String name_file = String.valueOf(gp.inter_dir) + gp.project_name + "_sample_names.txt";
        final DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
        System.out.println("PoolHapX " + gp.function + " Initiated: " + dtf.format(LocalDateTime.now()));
        final Evaluate eva = new Evaluate();
        if (function.equals("format")) {
            final String[] sam_files = get_filepaths(name_file, String.valueOf(gp.input_dir) + "/sam/", "sam", true);
            BAMFormatterGATK.generate_vef_varfreq(gp.input_dir, gp.inter_dir, gp.project_name, sam_files, (HashMap)Entrance.name2index, gp.sequencing_technology);
            final File file = new File(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_vars.intra_freq.txt");
            if (file.exists()) {
                BAMFormatterGATK.rewrite_inter_vars(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_vars.intra_freq.txt", gs_var_pos);
            }
        }
        else if (function.equals("gc")) {
            final String[] vef_files = get_filepaths(name_file, String.valueOf(gp.inter_dir) + "vef", "vef", false);
            new File(String.valueOf(gp.inter_dir) + "/gcf/").mkdir();
            for (int p = 0; p < Entrance.num_pools; ++p) {
                try {
                    final GraphColoring graphColoring = new GraphColoring(vef_files[p], gs_var_pos, String.valueOf(gp.inter_dir) + "gcf/" + Entrance.names_array[p] + ".gcf", gp.num_pos_window, gp.num_gap_window);
                }
                catch (Exception e3) {
                    final GraphColoring graphColoring2 = new GraphColoring(vef_files[p], gs_var_pos, String.valueOf(gp.inter_dir) + "gcf/" + Entrance.names_array[p] + ".gcf");
                }
                System.out.println("Graph colouring for pool " + p + ":" + Entrance.names_array[p] + " is finished.");
            }
            System.out.println("\nGraph-Coloring Finished: " + dtf.format(LocalDateTime.now()) + "\n");
        }
        else if (function.equals("aem_previous")) {
            final String dc_out_file = String.valueOf(gp.inter_dir) + gp.project_name + "_dc_plan.txt";
            final String[] vef_files2 = get_filepaths(name_file, String.valueOf(gp.inter_dir) + "vef", "vef", false);
            final String[] gcf_files = get_filepaths(name_file, String.valueOf(gp.inter_dir) + "gcf", "gcf", false);
            final DivideConquer dc_maker = new DivideConquer(gs_var_pos, gcf_files, parameter_file, dc_out_file);
            new File(String.valueOf(gp.inter_dir) + "/aem/").mkdir();
            final HapConfig[] level_I_config = dc_maker.regional_AEM(Entrance.names_array, vef_files2, dc_maker.regions_level_I, parameter_file, String.valueOf(gp.inter_dir) + "/aem/" + gp.project_name, String.valueOf(gp.inter_dir) + "/aem_fail_regional_lasso/", 1);
            System.out.println("Level 1 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
            final HapConfig[] level_II_config = dc_maker.regional_AEM(Entrance.names_array, vef_files2, dc_maker.regions_level_II, parameter_file, String.valueOf(gp.inter_dir) + "/aem/" + gp.project_name, String.valueOf(gp.inter_dir) + "/aem_fail_regional_lasso/", 2);
            System.out.println("Level 2 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
            eva.AemEvaluate(gp.project_name, dc_out_file, "/home/chencao/Desktop/sim001/gold_standard/sim001_haps.txt", "/home/chencao/Desktop/sim001/intermediate/aem/");
            final GraphColoring region_linker = new GraphColoring(level_I_config, level_II_config, gs_var_pos, gp.virtual_cov_link_gc);
            final HapConfig final_global_haps = region_linker.hapOut(Entrance.names_array);
            final_global_haps.recode_HapIDs_to_base16();
            final_global_haps.write_global_file_string(String.valueOf(gp.out_dir) + gp.project_name + "_gc.inter_freq_haps.txt");
            System.out.println("\nPost-AEM GC Finished: " + dtf.format(LocalDateTime.now()) + "\n");
            eva.GcAemEvaluate("/home/chencao/Desktop/sim001/gold_standard/sim001_haps.txt", "/home/chencao/Desktop/sim001/output/sim001_gc.inter_freq_haps.txt");
        }
        else if (function.equals("lasso")) {
            new File(String.valueOf(gp.inter_dir) + "/lasso/").mkdir();
            get_filepaths(name_file, String.valueOf(gp.inter_dir) + "/gcf/", "gcf", false);
            final HapConfig final_global_haps2 = new HapConfig(String.valueOf(gp.out_dir) + gp.project_name + "_hc.inter_freq_haps.txt", (String)null);
            final SiteInPoolFreqAnno siteInPoolFreqAnno = new SiteInPoolFreqAnno(gs_var_pos);
            final_global_haps2.inpool_site_freqs = siteInPoolFreqAnno.inpool_freqs;
            final_global_haps2.num_pools = Entrance.num_pools;
            final_global_haps2.pool_IDs = siteInPoolFreqAnno.pool_IDs;
            final_global_haps2.in_pool_haps_freq = new double[final_global_haps2.num_global_hap][final_global_haps2.num_pools];
            final HapConfig[] final_inpool_haps = new HapConfig[Entrance.num_pools];
            for (int pool_index = 0; pool_index < Entrance.num_pools; ++pool_index) {
                final HapLASSO inpool_lasso = new HapLASSO(pool_index, Entrance.names_array[pool_index], gp.lasso_global_lambda, final_global_haps2, gp.lasso_full_hap_freq_cutoff, gp.lasso_global_memory, String.valueOf(gp.inter_dir) + "/lasso/" + Entrance.names_array[pool_index], gp.lasso_coverage_weight, gp.lasso_distance_max_weight);
                inpool_lasso.estimate_frequencies_lasso(String.valueOf(gp.inter_dir) + "/vef/" + Entrance.names_array[pool_index] + ".vef", (String[])null, gp.lasso_weights);
            }
            System.out.println("\nGlobal LASSO Finished: " + dtf.format(LocalDateTime.now()) + "\n");
        }
        else if (function.equals("clustering")) {
            final HapConfig hapConfig = new HapConfig(String.valueOf(gp.out_dir) + gp.project_name + "_gc.inter_freq_haps.txt", (String)null);
        }
        else if (function.equals("aem")) {
            final String dc_out_file = String.valueOf(gp.inter_dir) + gp.project_name + "_dc_plan.txt";
            final String[] vef_files2 = get_filepaths(name_file, String.valueOf(gp.inter_dir) + "vef", "vef", false);
            final String[] gcf_files = get_filepaths(name_file, String.valueOf(gp.inter_dir) + "gcf", "gcf", false);
            final DivideConquer dc_maker = new DivideConquer(gs_var_pos, gcf_files, parameter_file, dc_out_file);
            new File(String.valueOf(gp.inter_dir) + "/aem/").mkdir();
            final HashMap<Integer, String> index_var_prefix_dict = (HashMap<Integer, String>)dc_maker.gs_map(gs_var_pos);
            final String reg_dc_out_file = String.valueOf(gp.inter_dir) + gp.project_name + "_regression__dc_plan.txt";
            final HashMap<Integer, Integer> index_var_pos_dict = (HashMap<Integer, Integer>)dc_maker.gs_map_pos(gs_var_pos);
            final DivideConquer reg_dc_maker = new DivideConquer(dc_maker.regions_level_final, dc_maker.regions_level_final_link, reg_dc_out_file, gp.regression_maximum_regions, (HashMap)index_var_pos_dict);
            System.out.println("DC Finished Time:\t" + dtf.format(LocalDateTime.now()) + "\n");
            for (int p2 = 0; p2 < Entrance.num_pools; ++p2) {
                final GraphColoring pool_in = new GraphColoring(vef_files2[p2], gs_var_pos, String.valueOf(gp.inter_dir) + "gcf/" + Entrance.names_array[p2] + ".gcf", dc_maker.regions_level_I);
                dc_maker.loci_link_freq[p2] = pool_in.loci_link_freq;
                System.out.println("Linking information extraction for pool " + p2 + ":" + Entrance.names_array[p2] + " is finished.");
            }
            HapConfig[] level_I_config2 = null;
            HapConfig[] level_II_config2 = null;
            HapConfig[] level_III_config = null;
            HapConfig[] level_IV_config = null;
            HapConfig[] level_V_config = null;
            HapConfig[] level_VI_config = null;
            HapConfig[] level_VII_config = null;
            HapConfig[] level_VIII_config = null;
            if (dc_maker.final_level >= 1) {
                level_I_config2 = dc_maker.regional_AEM(Entrance.names_array, vef_files2, dc_maker.regions_level_I, parameter_file, String.valueOf(gp.inter_dir) + "/aem/" + gp.project_name, String.valueOf(gp.inter_dir) + "/aem_fail_regional_lasso/", 1, (HashMap)Entrance.name2index);
                dc_maker.level_I_config = level_I_config2;
                System.out.println("Level 1 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
                level_II_config2 = dc_maker.regional_AEM(Entrance.names_array, vef_files2, dc_maker.regions_level_II, parameter_file, String.valueOf(gp.inter_dir) + "/aem/" + gp.project_name, String.valueOf(gp.inter_dir) + "/aem_fail_regional_lasso/", 2, (HashMap)Entrance.name2index);
                dc_maker.level_II_config = level_II_config2;
                System.out.println("Level 2 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
                eva.AemEvaluate(gp.project_name, dc_out_file, String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt", String.valueOf(gp.inter_dir) + "/aem/");
            }
            if (dc_maker.final_level >= 3) {
                level_III_config = dc_maker.level_III_regional_AEM(Entrance.names_array, vef_files2, dc_maker.regions_level_III, parameter_file, String.valueOf(gp.inter_dir) + "/aem/" + gp.project_name, String.valueOf(gp.inter_dir) + "/aem_fail_regional_lasso/", 3, (HashMap)Entrance.name2index, gp.level_III_IV_region_mismatch_tolerance);
                dc_maker.level_III_config = level_III_config;
                System.out.println("Level 3 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
                level_IV_config = dc_maker.level_IV_regional_AEM(Entrance.names_array, vef_files2, dc_maker.regions_level_IV, parameter_file, String.valueOf(gp.inter_dir) + "/aem/" + gp.project_name, String.valueOf(gp.inter_dir) + "/aem_fail_regional_lasso/", 4, (HashMap)Entrance.name2index, gp.level_III_IV_region_mismatch_tolerance);
                dc_maker.level_IV_config = level_IV_config;
                System.out.println("Level 4 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
                eva.Aem_III_IV_Evaluate(gp.project_name, dc_out_file, String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt", String.valueOf(gp.inter_dir) + "/aem/");
            }
            if (dc_maker.final_level >= 5) {
                level_V_config = dc_maker.level_V_regional_AEM(Entrance.names_array, vef_files2, dc_maker.regions_level_V, parameter_file, String.valueOf(gp.inter_dir) + "/aem/" + gp.project_name, String.valueOf(gp.inter_dir) + "/aem_fail_regional_lasso/", 5, (HashMap)Entrance.name2index, gp.level_V_VI_region_mismatch_tolerance);
                dc_maker.level_V_config = level_V_config;
                System.out.println("Level 5 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
                level_VI_config = dc_maker.level_VI_regional_AEM(Entrance.names_array, vef_files2, dc_maker.regions_level_VI, parameter_file, String.valueOf(gp.inter_dir) + "/aem/" + gp.project_name, String.valueOf(gp.inter_dir) + "/aem_fail_regional_lasso/", 6, (HashMap)Entrance.name2index, gp.level_V_VI_region_mismatch_tolerance);
                dc_maker.level_VI_config = level_VI_config;
                System.out.println("Level 6 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
                eva.Aem_V_VI_Evaluate(gp.project_name, dc_out_file, String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt", String.valueOf(gp.inter_dir) + "/aem/");
            }
            if (dc_maker.final_level >= 7) {
                level_VII_config = dc_maker.level_VII_regional_AEM(Entrance.names_array, vef_files2, dc_maker.regions_level_VII, parameter_file, String.valueOf(gp.inter_dir) + "/aem/" + gp.project_name, String.valueOf(gp.inter_dir) + "/aem_fail_regional_lasso/", 7, (HashMap)Entrance.name2index, gp.level_VII_VIII_region_mismatch_tolerance);
                dc_maker.level_VII_config = level_VII_config;
                System.out.println("Level 7 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
                level_VIII_config = dc_maker.level_VIII_regional_AEM(Entrance.names_array, vef_files2, dc_maker.regions_level_VIII, parameter_file, String.valueOf(gp.inter_dir) + "/aem/" + gp.project_name, String.valueOf(gp.inter_dir) + "/aem_fail_regional_lasso/", 8, (HashMap)Entrance.name2index, gp.level_VII_VIII_region_mismatch_tolerance);
                dc_maker.level_VIII_config = level_VIII_config;
                System.out.println("Level 8 AEM Finished: " + dtf.format(LocalDateTime.now()) + "\n");
                eva.Aem_VII_VIII_Evaluate(gp.project_name, dc_out_file, String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt", String.valueOf(gp.inter_dir) + "/aem/");
            }
            dc_maker.level_final_config = null;
            dc_maker.level_final_link_config = null;
            if (dc_maker.final_level == 7) {
                dc_maker.level_final_config = level_VII_config;
                dc_maker.level_final_link_config = level_VIII_config;
            }
            else if (dc_maker.final_level == 5) {
                dc_maker.level_final_config = level_V_config;
                dc_maker.level_final_link_config = level_VI_config;
            }
            else if (dc_maker.final_level == 3) {
                dc_maker.level_final_config = level_III_config;
                dc_maker.level_final_link_config = level_IV_config;
            }
            else if (dc_maker.final_level == 1) {
                dc_maker.level_final_config = level_I_config2;
                dc_maker.level_final_link_config = level_II_config2;
            }
            System.out.println("AEM Finished Time:\t" + dtf.format(LocalDateTime.now()) + "\n");
            System.out.println("\nAEM Finished and Start regional L0L1 regression: " + dtf.format(LocalDateTime.now()) + "\n");
            final String reg_prefix = String.valueOf(gp.inter_dir) + "/regression/";
            try {
                final String shpath = "rm -rf  " + reg_prefix;
                final Process ps = Runtime.getRuntime().exec(shpath, null);
                ps.waitFor();
                final BufferedReader br = new BufferedReader(new InputStreamReader(ps.getInputStream()));
                final StringBuffer sb = new StringBuffer();
                String line;
                while ((line = br.readLine()) != null) {
                    sb.append(line).append("\n");
                }
                final String result = sb.toString();
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            new File(reg_prefix).mkdir();
            new File(String.valueOf(reg_prefix) + "/Rscript").mkdir();
            for (int pool_index2 = 0; pool_index2 < Entrance.num_pools; ++pool_index2) {
                new File(String.valueOf(gp.inter_dir) + "/regression/" + Entrance.names_array[pool_index2]).mkdir();
                for (int r = 0; r < reg_dc_maker.regession_level[0].length; ++r) {
                    final GraphColoring graphColoring3 = new GraphColoring(reg_dc_maker.regession_level[0][r][0], reg_dc_maker.regession_level[0][r][1], dc_maker.level_final_config, dc_maker.level_final_link_config, gs_var_pos, dc_maker.regions_level_final, dc_maker.regions_level_final_link, gp.bfs_mismatch_tolerance, r + 1, (HashMap)index_var_prefix_dict, String.valueOf(reg_prefix) + Entrance.names_array[pool_index2], gp.regression_maximum_regions);
                }
            }
            System.out.println("Have generated initialization regional haptypes from AEM level " + dc_maker.final_level + " regional haplotypes. " + dtf.format(LocalDateTime.now()) + "\n");
            for (int regression_level = 0; regression_level < reg_dc_maker.regession_level.length; ++regression_level) {
                for (int r = 0; r < reg_dc_maker.regession_level[regression_level].length; ++r) {
                    for (int pool_index3 = 0; pool_index3 < Entrance.num_pools; ++pool_index3) {
                        System.out.println("Generating level " + Integer.toString(regression_level + 1) + " region " + Integer.toString(r + 1) + "  potential haplotypes for pool " + Entrance.names_array[pool_index3] + ". " + dtf.format(LocalDateTime.now()));
                        if (reg_dc_maker.regession_level[regression_level][r][0] != -1 && reg_dc_maker.regession_level[regression_level][r][1] != -1) {
                            if (regression_level != 0) {
                                final GraphColoring graphColoring4 = new GraphColoring(reg_dc_maker.regession_level[regression_level][r][0], reg_dc_maker.regession_level[regression_level][r][1], reg_dc_maker.regession_level[regression_level - 1], dc_maker.level_final_link_config, gs_var_pos, dc_maker.regions_level_final_link, gp.regression_link_mismatch_tolerance, r + 1, regression_level + 1, (HashMap)index_var_prefix_dict, String.valueOf(reg_prefix) + Entrance.names_array[pool_index3], gp.regression_maximum_regions, gp.regression_maximum_selected_haplotypes);
                            }
                            final HapConfig final_global_haps3 = new HapConfig(String.valueOf(gp.inter_dir) + "/regression/" + Entrance.names_array[pool_index3] + "/regression_level_" + Integer.toString(regression_level + 1) + "_region_" + Integer.toString(r + 1) + ".potential.haps", (String)null);
                            final SiteInPoolFreqAnno siteInPoolFreqAnno2 = new SiteInPoolFreqAnno(gs_var_pos, reg_dc_maker.regession_level[regression_level][r][0], reg_dc_maker.regession_level[regression_level][r][1], pool_index3);
                            final_global_haps3.inpool_site_freqs = siteInPoolFreqAnno2.inpool_freqs;
                            final_global_haps3.num_pools = Entrance.num_pools;
                            final_global_haps3.pool_IDs = siteInPoolFreqAnno2.pool_IDs;
                            final_global_haps3.in_pool_haps_freq = new double[final_global_haps3.num_global_hap][final_global_haps3.num_pools];
                            final HapLASSO inpool_lasso2 = new HapLASSO(pool_index3, Entrance.names_array[pool_index3], gp.lasso_global_lambda, final_global_haps3, gp.lasso_full_hap_freq_cutoff, gp.lasso_global_memory, String.valueOf(gp.inter_dir) + "/regression/" + Entrance.names_array[pool_index3] + "/regression_level_" + Integer.toString(regression_level + 1) + "_region_" + Integer.toString(r + 1), gp.lasso_coverage_weight, gp.lasso_distance_max_weight);
                            inpool_lasso2.estimate_frequencies_lasso(String.valueOf(gp.inter_dir) + "/vef/" + Entrance.names_array[pool_index3] + ".vef", (String[])null, gp.lasso_weights, gp.sequencing_technology, (HashMap)reg_dc_maker.regession_level_pos_come_from.get(regression_level).get(r), reg_dc_maker.regession_level_num_come_from[regression_level][r], (HashMap)index_var_pos_dict, regression_level);
                        }
                    }
                }
                System.out.println("\nRunning L0L1 regression for all regions in level " + Integer.toString(regression_level + 1) + ". " + dtf.format(LocalDateTime.now()) + "\n");
                final HapLASSO regression_rscript = new HapLASSO(Entrance.names_array, reg_dc_maker.regession_level[regression_level], String.valueOf(gp.inter_dir) + "/regression/Rscript/", String.valueOf(gp.inter_dir) + "/regression/", regression_level + 1, gp.regression_hapset_size_max, gp.regression_hapset_size_min, gp.regression_gamma_min, gp.regression_gamma_max, gp.ngammas, gp.lasso_weights[1], gp.species);
                for (int r2 = 0; r2 < reg_dc_maker.regession_level[regression_level].length; ++r2) {
                    if (reg_dc_maker.regession_level[regression_level][r2][0] != -1 && reg_dc_maker.regession_level[regression_level][r2][1] != -1) {
                        regression_rscript.run_L0L1learn(Entrance.names_array, gp.rscript_path, String.valueOf(gp.inter_dir) + "/regression/Rscript/regression_level_" + Integer.toString(regression_level + 1) + "_region_" + Integer.toString(r2 + 1), String.valueOf(gp.inter_dir) + "/regression/", regression_level + 1, r2 + 1, gp.num_threads);
                        System.out.println("Level " + Integer.toString(regression_level + 1) + " region " + Integer.toString(r2 + 1) + " L0L1 regression finished.\t" + dtf.format(LocalDateTime.now()) + "\n");
                    }
                }
            }
            final int regression_last_level = reg_dc_maker.regession_level.length - 1;
            for (int pool_index4 = 0; pool_index4 < Entrance.num_pools; ++pool_index4) {
                new File(String.valueOf(gp.out_dir) + Entrance.names_array[pool_index4]).mkdir();
                final HapConfig hapConfig2 = new HapConfig(String.valueOf(gp.inter_dir) + "/regression/" + Entrance.names_array[pool_index4] + "/regression_level_" + Integer.toString(regression_last_level + 1) + "_region_1.regression_out", String.valueOf(gp.out_dir) + Entrance.names_array[pool_index4] + "/final_freq_haps.txt", (HashMap)index_var_prefix_dict);
            }
            System.out.println("L0L1 Regression Finished Time:\t" + dtf.format(LocalDateTime.now()) + "\n");
            System.out.println("PoolHapX Successfully Finished, Enjoy!\n");
        }
        else if (function.equals("evaluate")) {
            final String[] vef_files = get_filepaths(name_file, String.valueOf(gp.inter_dir) + "vef", "vef", false);
            final BufferedWriter bw_mcc = new BufferedWriter(new FileWriter(String.valueOf(gp.out_dir) + "/MCC.result"));
            double mcc_total = 0.0;
            for (int pool_index5 = 0; pool_index5 < Entrance.num_pools; ++pool_index5) {
                eva.MCCEvaluate(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt", String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.intra_freq.txt", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index5] + "/final_freq_haps.txt", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index5] + "/MCC.txt", gp.mcc_freq_cutoff, Entrance.names_array[pool_index5]);
                bw_mcc.write("MCC for " + Entrance.names_array[pool_index5] + " is:\t" + eva.mcc_value + "\n");
                mcc_total += eva.mcc_value;
            }
            System.out.println("Average MCC for all " + Entrance.num_pools + " pools:\t" + mcc_total / Entrance.num_pools);
            bw_mcc.write("Average MCC for all " + Entrance.num_pools + " pools:\t" + mcc_total / Entrance.num_pools + "\n");
            bw_mcc.close();
            final BufferedWriter bw_jsd = new BufferedWriter(new FileWriter(String.valueOf(gp.out_dir) + "/JSD.result"));
            double jsd_total = 0.0;
            for (int pool_index6 = 0; pool_index6 < Entrance.num_pools; ++pool_index6) {
                eva.JSDEvaluate(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt", String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.intra_freq.txt", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index6] + "/final_freq_haps.txt", Entrance.names_array[pool_index6], String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index6] + "/JSD.txt");
                jsd_total += eva.jsd_value;
                bw_jsd.write("JSD for " + Entrance.names_array[pool_index6] + " is:\t" + eva.jsd_value + "\n");
            }
            System.out.println("Average JSD for all " + Entrance.num_pools + " pools:\t" + jsd_total / Entrance.num_pools);
            bw_jsd.write("Average JSD for all " + Entrance.num_pools + " pools:\t" + jsd_total / Entrance.num_pools + "\n");
            bw_jsd.close();
            System.out.println("\nEvaluation Finished: " + dtf.format(LocalDateTime.now()) + "\n");
        }
        else if (function.equals("l0l1")) {
            final String dc_out_file = String.valueOf(gp.inter_dir) + gp.project_name + "_dc_plan.txt";
            final String[] vef_files2 = get_filepaths(name_file, String.valueOf(gp.inter_dir) + "vef", "vef", false);
            final String[] gcf_files = get_filepaths(name_file, String.valueOf(gp.inter_dir) + "gcf", "gcf", false);
            final DivideConquer dc_maker = new DivideConquer(gs_var_pos, gcf_files, parameter_file, dc_out_file);
            new File(String.valueOf(gp.inter_dir) + "/aem/").mkdir();
            final HashMap<Integer, String> index_var_prefix_dict = (HashMap<Integer, String>)dc_maker.gs_map(gs_var_pos);
            final String reg_dc_out_file = String.valueOf(gp.inter_dir) + gp.project_name + "_regression__dc_plan.txt";
            final HashMap<Integer, Integer> index_var_pos_dict = (HashMap<Integer, Integer>)dc_maker.gs_map_pos(gs_var_pos);
            final DivideConquer reg_dc_maker = new DivideConquer(dc_maker.regions_level_final, dc_maker.regions_level_final_link, reg_dc_out_file, gp.regression_maximum_regions, (HashMap)index_var_pos_dict);
            System.out.println("DC Finished:\t" + dtf.format(LocalDateTime.now()) + "\n");
            dc_maker.level_final_config = new HapConfig[dc_maker.num_regions_level_final];
            dc_maker.level_final_link_config = new HapConfig[dc_maker.num_regions_level_final_link];
            for (int r3 = 0; r3 < dc_maker.num_regions_level_final; ++r3) {
                final HapConfig final_aem_config = new HapConfig(String.valueOf(gp.inter_dir) + "aem/" + gp.project_name + "_level_" + Integer.toString(dc_maker.final_level) + "_region_" + r3 + ".inter_freq_haps.txt", (String)null);
                dc_maker.level_final_config[r3] = final_aem_config;
            }
            for (int r3 = 0; r3 < dc_maker.num_regions_level_final_link; ++r3) {
                final HapConfig final_link_aem_config = new HapConfig(String.valueOf(gp.inter_dir) + "aem/" + gp.project_name + "_level_" + Integer.toString(dc_maker.final_level + 1) + "_region_" + r3 + ".inter_freq_haps.txt", (String)null);
                dc_maker.level_final_link_config[r3] = final_link_aem_config;
            }
            System.out.println("\nAEM Finished and Start regional L0L1 regression: " + dtf.format(LocalDateTime.now()) + "\n");
            final String reg_prefix2 = String.valueOf(gp.inter_dir) + "/regression/";
            try {
                final String shpath2 = "rm -rf  " + reg_prefix2;
                final Process ps2 = Runtime.getRuntime().exec(shpath2, null);
                ps2.waitFor();
                final BufferedReader br2 = new BufferedReader(new InputStreamReader(ps2.getInputStream()));
                final StringBuffer sb2 = new StringBuffer();
                String line2;
                while ((line2 = br2.readLine()) != null) {
                    sb2.append(line2).append("\n");
                }
                final String result2 = sb2.toString();
                br2.close();
            }
            catch (Exception e2) {
                e2.printStackTrace();
            }
            new File(reg_prefix2).mkdir();
            new File(String.valueOf(reg_prefix2) + "/Rscript").mkdir();
            for (int pool_index7 = 0; pool_index7 < Entrance.num_pools; ++pool_index7) {
                new File(String.valueOf(gp.inter_dir) + "/regression/" + Entrance.names_array[pool_index7]).mkdir();
                for (int r4 = 0; r4 < reg_dc_maker.regession_level[0].length; ++r4) {
                    final GraphColoring graphColoring5 = new GraphColoring(reg_dc_maker.regession_level[0][r4][0], reg_dc_maker.regession_level[0][r4][1], dc_maker.level_final_config, dc_maker.level_final_link_config, gs_var_pos, dc_maker.regions_level_final, dc_maker.regions_level_final_link, gp.bfs_mismatch_tolerance, r4 + 1, (HashMap)index_var_prefix_dict, String.valueOf(reg_prefix2) + Entrance.names_array[pool_index7], gp.regression_maximum_regions);
                }
            }
            System.out.println("Have generated initialization regional haptypes from AEM level " + dc_maker.final_level + " regional haplotypes. " + dtf.format(LocalDateTime.now()) + "\n");
            for (int regression_level2 = 0; regression_level2 < reg_dc_maker.regession_level.length; ++regression_level2) {
                for (int r4 = 0; r4 < reg_dc_maker.regession_level[regression_level2].length; ++r4) {
                    for (int pool_index8 = 0; pool_index8 < Entrance.num_pools; ++pool_index8) {
                        System.out.println("Generating level " + Integer.toString(regression_level2 + 1) + " region " + Integer.toString(r4 + 1) + "  potential haplotypes for pool " + Entrance.names_array[pool_index8] + ". " + dtf.format(LocalDateTime.now()));
                        if (reg_dc_maker.regession_level[regression_level2][r4][0] != -1 && reg_dc_maker.regession_level[regression_level2][r4][1] != -1) {
                            if (regression_level2 != 0) {
                                final GraphColoring graphColoring6 = new GraphColoring(reg_dc_maker.regession_level[regression_level2][r4][0], reg_dc_maker.regession_level[regression_level2][r4][1], reg_dc_maker.regession_level[regression_level2 - 1], dc_maker.level_final_link_config, gs_var_pos, dc_maker.regions_level_final_link, gp.regression_link_mismatch_tolerance, r4 + 1, regression_level2 + 1, (HashMap)index_var_prefix_dict, String.valueOf(reg_prefix2) + Entrance.names_array[pool_index8], gp.regression_maximum_regions, gp.regression_maximum_selected_haplotypes);
                            }
                            final HapConfig final_global_haps4 = new HapConfig(String.valueOf(gp.inter_dir) + "/regression/" + Entrance.names_array[pool_index8] + "/regression_level_" + Integer.toString(regression_level2 + 1) + "_region_" + Integer.toString(r4 + 1) + ".potential.haps", (String)null);
                            final SiteInPoolFreqAnno siteInPoolFreqAnno3 = new SiteInPoolFreqAnno(gs_var_pos, reg_dc_maker.regession_level[regression_level2][r4][0], reg_dc_maker.regession_level[regression_level2][r4][1], pool_index8);
                            final_global_haps4.inpool_site_freqs = siteInPoolFreqAnno3.inpool_freqs;
                            final_global_haps4.num_pools = Entrance.num_pools;
                            final_global_haps4.pool_IDs = siteInPoolFreqAnno3.pool_IDs;
                            final_global_haps4.in_pool_haps_freq = new double[final_global_haps4.num_global_hap][final_global_haps4.num_pools];
                            final HapLASSO inpool_lasso3 = new HapLASSO(pool_index8, Entrance.names_array[pool_index8], gp.lasso_global_lambda, final_global_haps4, gp.lasso_full_hap_freq_cutoff, gp.lasso_global_memory, String.valueOf(gp.inter_dir) + "/regression/" + Entrance.names_array[pool_index8] + "/regression_level_" + Integer.toString(regression_level2 + 1) + "_region_" + Integer.toString(r4 + 1), gp.lasso_coverage_weight, gp.lasso_distance_max_weight);
                            inpool_lasso3.estimate_frequencies_lasso(String.valueOf(gp.inter_dir) + "/vef/" + Entrance.names_array[pool_index8] + ".vef", (String[])null, gp.lasso_weights, gp.sequencing_technology, (HashMap)reg_dc_maker.regession_level_pos_come_from.get(regression_level2).get(r4), reg_dc_maker.regession_level_num_come_from[regression_level2][r4], (HashMap)index_var_pos_dict, regression_level2);
                        }
                    }
                }
                System.out.println("\nRunning L0L1 regression for all regions in level " + Integer.toString(regression_level2 + 1) + ". " + dtf.format(LocalDateTime.now()) + "\n");
                final HapLASSO regression_rscript2 = new HapLASSO(Entrance.names_array, reg_dc_maker.regession_level[regression_level2], String.valueOf(gp.inter_dir) + "/regression/Rscript/", String.valueOf(gp.inter_dir) + "/regression/", regression_level2 + 1, gp.regression_hapset_size_max, gp.regression_hapset_size_min, gp.regression_gamma_min, gp.regression_gamma_max, gp.ngammas, gp.lasso_weights[1], gp.species);
                for (int r5 = 0; r5 < reg_dc_maker.regession_level[regression_level2].length; ++r5) {
                    if (reg_dc_maker.regession_level[regression_level2][r5][0] != -1 && reg_dc_maker.regession_level[regression_level2][r5][1] != -1) {
                        regression_rscript2.run_L0L1learn(Entrance.names_array, gp.rscript_path, String.valueOf(gp.inter_dir) + "/regression/Rscript/regression_level_" + Integer.toString(regression_level2 + 1) + "_region_" + Integer.toString(r5 + 1), String.valueOf(gp.inter_dir) + "/regression/", regression_level2 + 1, r5 + 1, gp.num_threads);
                        System.out.println("Level " + Integer.toString(regression_level2 + 1) + " region " + Integer.toString(r5 + 1) + " L0L1 regression finished.\t" + dtf.format(LocalDateTime.now()) + "\n");
                    }
                }
            }
            final int regression_last_level2 = reg_dc_maker.regession_level.length - 1;
            for (int pool_index9 = 0; pool_index9 < Entrance.num_pools; ++pool_index9) {
                new File(String.valueOf(gp.out_dir) + Entrance.names_array[pool_index9]).mkdir();
                final HapConfig hapConfig3 = new HapConfig(String.valueOf(gp.inter_dir) + "/regression/" + Entrance.names_array[pool_index9] + "/regression_level_" + Integer.toString(regression_last_level2 + 1) + "_region_1.regression_out", String.valueOf(gp.out_dir) + Entrance.names_array[pool_index9] + "/final_freq_haps.txt", (HashMap)index_var_prefix_dict);
            }
            System.out.println("PoolHapX Successfully Finished, Enjoy!\n");
        }
        else if (function.equals("Comparison_bacteria")) {
            final String[] vef_files = get_filepaths(name_file, String.valueOf(gp.inter_dir) + "vef", "vef", false);
            final String dc_out_file2 = String.valueOf(gp.inter_dir) + gp.project_name + "_dc_plan.txt";
            final String[] gcf_files = get_filepaths(name_file, String.valueOf(gp.inter_dir) + "gcf", "gcf", false);
            final DivideConquer dc_maker = new DivideConquer();
            final HashMap<Integer, String> index_var_prefix_dict = (HashMap<Integer, String>)dc_maker.gs_map(gs_var_pos);
            final BufferedWriter bw_bhap_mcc = new BufferedWriter(new FileWriter(String.valueOf(gp.out_dir) + "/MCC_bhap.result"));
            double mcc_total2 = 0.0;
            for (int pool_index10 = 0; pool_index10 < Entrance.num_pools; ++pool_index10) {
                final String[] tmp = Entrance.names_array[pool_index10].split("p");
                final String pool_index_name = tmp[tmp.length - 1];
                eva.GenerateFinal_BHap(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_vars.intra_freq.txt", String.valueOf(gp.out_dir) + "/bhap_pool_" + pool_index_name + "/finalPredictions", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index10] + "/bhap_final_freq_haps.txt", (HashMap)index_var_prefix_dict, Entrance.names_array[pool_index10]);
                eva.MCCEvaluate(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt", String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.intra_freq.txt", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index10] + "/bhap_final_freq_haps.txt", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index10] + "/bhap_MCC.txt", gp.mcc_freq_cutoff, Entrance.names_array[pool_index10]);
                bw_bhap_mcc.write("MCC for " + Entrance.names_array[pool_index10] + " is:\t" + eva.mcc_value + "\n");
                mcc_total2 += eva.mcc_value;
            }
            System.out.println("Average bhap MCC for all " + Entrance.num_pools + " pools:\t" + mcc_total2 / Entrance.num_pools);
            bw_bhap_mcc.write("Average bhap MCC for all " + Entrance.num_pools + " pools:\t" + mcc_total2 / Entrance.num_pools + "\n");
            bw_bhap_mcc.close();
            final BufferedWriter bw_bhap_jsd = new BufferedWriter(new FileWriter(String.valueOf(gp.out_dir) + "/bhap_JSD.result"));
            double jsd_total2 = 0.0;
            for (int pool_index8 = 0; pool_index8 < Entrance.num_pools; ++pool_index8) {
                eva.JSDEvaluate(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt", String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.intra_freq.txt", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index8] + "/bhap_final_freq_haps.txt", Entrance.names_array[pool_index8], String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index8] + "/bhap_JSD.txt");
                jsd_total2 += eva.jsd_value;
                bw_bhap_jsd.write("JSD for " + Entrance.names_array[pool_index8] + " is:\t" + eva.jsd_value + "\n");
            }
            System.out.println("Average bhap JSD for all " + Entrance.num_pools + " pools:\t" + jsd_total2 / Entrance.num_pools);
            bw_bhap_jsd.write("Average bhap JSD for all " + Entrance.num_pools + " pools:\t" + jsd_total2 / Entrance.num_pools + "\n");
            bw_bhap_jsd.close();
            final BufferedWriter bw_evorha_mcc = new BufferedWriter(new FileWriter(String.valueOf(gp.out_dir) + "/MCC_evorha.result"));
            mcc_total2 = 0.0;
            for (int pool_index11 = 0; pool_index11 < Entrance.num_pools; ++pool_index11) {
                final String[] tmp2 = Entrance.names_array[pool_index11].split("p");
                final String pool_index_name2 = tmp2[tmp2.length - 1];
                eva.GenerateFinal_EVORHA(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_mutations.txt", String.valueOf(gp.out_dir) + "/evorha_pool_" + pool_index_name2 + "/evorha.global.hapfreq", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index11] + "/evorha_final_freq_haps.txt", (HashMap)index_var_prefix_dict, Entrance.names_array[pool_index11]);
                eva.MCCEvaluate(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt", String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.intra_freq.txt", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index11] + "/evorha_final_freq_haps.txt", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index11] + "/evorha_MCC.txt", gp.mcc_freq_cutoff, Entrance.names_array[pool_index11]);
                bw_evorha_mcc.write("MCC for " + Entrance.names_array[pool_index11] + " is:\t" + eva.mcc_value + "\n");
                mcc_total2 += eva.mcc_value;
            }
            System.out.println("Average evorha MCC for all " + Entrance.num_pools + " pools:\t" + mcc_total2 / Entrance.num_pools);
            bw_evorha_mcc.write("Average evorha MCC for all " + Entrance.num_pools + " pools:\t" + mcc_total2 / Entrance.num_pools + "\n");
            bw_evorha_mcc.close();
            final BufferedWriter bw_evorha_jsd = new BufferedWriter(new FileWriter(String.valueOf(gp.out_dir) + "/evorha_JSD.result"));
            jsd_total2 = 0.0;
            for (int pool_index12 = 0; pool_index12 < Entrance.num_pools; ++pool_index12) {
                eva.JSDEvaluate(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt", String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.intra_freq.txt", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index12] + "/evorha_final_freq_haps.txt", Entrance.names_array[pool_index12], String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index12] + "/evorha_JSD.txt");
                jsd_total2 += eva.jsd_value;
                bw_evorha_jsd.write("JSD for " + Entrance.names_array[pool_index12] + " is:\t" + eva.jsd_value + "\n");
            }
            System.out.println("Average evorha JSD for all " + Entrance.num_pools + " pools:\t" + jsd_total2 / Entrance.num_pools);
            bw_evorha_jsd.write("Average evorha JSD for all " + Entrance.num_pools + " pools:\t" + jsd_total2 / Entrance.num_pools + "\n");
            bw_evorha_jsd.close();
        }
        else if (function.equals("Comparison_metagenomics")) {
            final String[] vef_files = get_filepaths(name_file, String.valueOf(gp.inter_dir) + "vef", "vef", false);
            final String dc_out_file2 = String.valueOf(gp.inter_dir) + gp.project_name + "_dc_plan.txt";
            final String[] gcf_files = get_filepaths(name_file, String.valueOf(gp.inter_dir) + "gcf", "gcf", false);
            final DivideConquer dc_maker = new DivideConquer();
            final HashMap<Integer, String> index_var_prefix_dict = (HashMap<Integer, String>)dc_maker.gs_map(gs_var_pos);
            final BufferedWriter bw_strainest_mcc = new BufferedWriter(new FileWriter(String.valueOf(gp.out_dir) + "/MCC_strainest.result"));
            double mcc_total2 = 0.0;
            for (int pool_index10 = 0; pool_index10 < Entrance.num_pools; ++pool_index10) {
                final String[] tmp = Entrance.names_array[pool_index10].split("p");
                final String pool_index_name = tmp[tmp.length - 1];
                eva.GenerateFinal_Strainest(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_vars.intra_freq.txt", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index10] + "/strainest_final_freq_haps.txt", (HashMap)index_var_prefix_dict, Entrance.names_array[pool_index10]);
                eva.MCCEvaluate(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt", String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.intra_freq.txt", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index10] + "/strainest_final_freq_haps.txt", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index10] + "/strainest_MCC.txt", gp.mcc_freq_cutoff, Entrance.names_array[pool_index10]);
                bw_strainest_mcc.write("MCC for " + Entrance.names_array[pool_index10] + " is:\t" + eva.mcc_value + "\n");
                mcc_total2 += eva.mcc_value;
            }
            System.out.println("Average strainest MCC for all " + Entrance.num_pools + " pools:\t" + mcc_total2 / Entrance.num_pools);
            bw_strainest_mcc.write("Average strainest MCC for all " + Entrance.num_pools + " pools:\t" + mcc_total2 / Entrance.num_pools + "\n");
            bw_strainest_mcc.close();
            final BufferedWriter bw_strainest_jsd = new BufferedWriter(new FileWriter(String.valueOf(gp.out_dir) + "/strainest_JSD.result"));
            double jsd_total2 = 0.0;
            for (int pool_index8 = 0; pool_index8 < Entrance.num_pools; ++pool_index8) {
                eva.JSDEvaluate(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt", String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.intra_freq.txt", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index8] + "/strainest_final_freq_haps.txt", Entrance.names_array[pool_index8], String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index8] + "/strainest_JSD.txt");
                jsd_total2 += eva.jsd_value;
                bw_strainest_jsd.write("JSD for " + Entrance.names_array[pool_index8] + " is:\t" + eva.jsd_value + "\n");
            }
            System.out.println("Average strainest JSD for all " + Entrance.num_pools + " pools:\t" + jsd_total2 / Entrance.num_pools);
            bw_strainest_jsd.write("Average strainest JSD for all " + Entrance.num_pools + " pools:\t" + jsd_total2 / Entrance.num_pools + "\n");
            bw_strainest_jsd.close();
            final BufferedWriter bw_gretel_mcc = new BufferedWriter(new FileWriter(String.valueOf(gp.out_dir) + "/MCC_gretel.result"));
            mcc_total2 = 0.0;
            double real_num_pools = 1.0E-5;
            for (int pool_index13 = 0; pool_index13 < Entrance.num_pools; ++pool_index13) {
                final String[] tmp3 = Entrance.names_array[pool_index13].split("p");
                final String pool_index_name3 = tmp3[tmp3.length - 1];
                final File file2 = new File(String.valueOf(gp.out_dir) + "/gretel_pool_" + pool_index_name3 + "/snp.fasta");
                if (file2.exists()) {
                    ++real_num_pools;
                    eva.GenerateFinal_Gretel(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_mutations.txt", String.valueOf(gp.out_dir) + "/gretel_pool_" + pool_index_name3 + "/gretel.vcf", String.valueOf(gp.out_dir) + "/gretel_pool_" + pool_index_name3 + "/snp.fasta", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index13] + "/gretel_final_freq_haps.txt", (HashMap)index_var_prefix_dict, Entrance.names_array[pool_index13]);
                    eva.MCCEvaluate(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt", String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.intra_freq.txt", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index13] + "/gretel_final_freq_haps.txt", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index13] + "/gretel_MCC.txt", gp.mcc_freq_cutoff, Entrance.names_array[pool_index13]);
                    bw_gretel_mcc.write("MCC for " + Entrance.names_array[pool_index13] + " is:\t" + eva.mcc_value + "\n");
                    mcc_total2 += eva.mcc_value;
                }
            }
            System.out.println("Average gretel MCC for all " + Entrance.num_pools + " pools:\t" + mcc_total2 / real_num_pools);
            bw_gretel_mcc.write("Average gretel MCC for all " + Entrance.num_pools + " pools:\t" + mcc_total2 / real_num_pools + "\n");
            bw_gretel_mcc.close();
            final BufferedWriter bw_gretel_jsd = new BufferedWriter(new FileWriter(String.valueOf(gp.out_dir) + "/gretel_JSD.result"));
            jsd_total2 = 0.0;
            for (int pool_index14 = 0; pool_index14 < Entrance.num_pools; ++pool_index14) {
                final String[] tmp4 = Entrance.names_array[pool_index14].split("p");
                final String pool_index_name4 = tmp4[tmp4.length - 1];
                final File file3 = new File(String.valueOf(gp.out_dir) + "/gretel_pool_" + pool_index_name4 + "/snp.fasta");
                if (file3.exists()) {
                    eva.JSDEvaluate(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt", String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.intra_freq.txt", String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index14] + "/gretel_final_freq_haps.txt", Entrance.names_array[pool_index14], String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index14] + "/gretel_JSD.txt");
                    jsd_total2 += eva.jsd_value;
                    bw_gretel_jsd.write("JSD for " + Entrance.names_array[pool_index14] + " is:\t" + eva.jsd_value + "\n");
                }
            }
            System.out.println("Average gretel JSD for all " + Entrance.num_pools + " pools:\t" + jsd_total2 / real_num_pools);
            bw_gretel_jsd.write("Average gretel JSD for all " + Entrance.num_pools + " pools:\t" + jsd_total2 / real_num_pools + "\n");
            bw_gretel_jsd.close();
        }
        else if (function.equals("Comparison_human")) {
            int default_sel_haps = 50;
            final String[] vef_files2 = get_filepaths(name_file, String.valueOf(gp.inter_dir) + "vef", "vef", false);
            final String dc_out_file3 = String.valueOf(gp.inter_dir) + gp.project_name + "_dc_plan.txt";
            final String[] gcf_files2 = get_filepaths(name_file, String.valueOf(gp.inter_dir) + "gcf", "gcf", false);
            final DivideConquer dc_maker2 = new DivideConquer();
            final HashMap<Integer, String> index_var_prefix_dict2 = (HashMap<Integer, String>)dc_maker2.gs_map(gs_var_pos);
            double freq_cutoff = 0.0;
            if (index_var_prefix_dict2.size() > 8 && index_var_prefix_dict2.size() < 12) {
                freq_cutoff = 0.02;
            }
            if (index_var_prefix_dict2.size() > 13 && index_var_prefix_dict2.size() < 17) {
                freq_cutoff = 0.015;
            }
            if (index_var_prefix_dict2.size() > 23 && index_var_prefix_dict2.size() < 27) {
                freq_cutoff = 0.01;
            }
            if (index_var_prefix_dict2.size() > 100 && index_var_prefix_dict2.size() < 10000) {
                freq_cutoff = 0.0;
                default_sel_haps = 500;
            }
            final String freq_fil = String.valueOf(gp.gold_dir) + gp.project_name + "_haps.inter_freq_vars.txt";
            final String[] tmp_arr1 = freq_fil.split("/");
            final String folder_prefix = tmp_arr1[tmp_arr1.length - 6];
            final String[] tmp_arr2 = folder_prefix.split("_");
            final String folder = String.valueOf(tmp_arr2[0]) + "_10";
            String new_freq_fil = tmp_arr1[0];
            for (int i = 1; i < tmp_arr1.length; ++i) {
                if (i != tmp_arr1.length - 6) {
                    new_freq_fil = String.valueOf(new_freq_fil) + "/" + tmp_arr1[i];
                }
                else {
                    new_freq_fil = String.valueOf(new_freq_fil) + "/" + folder;
                }
            }
            System.out.println(new_freq_fil);
            final String godl_freq_file = new_freq_fil;
            final BufferedReader br_g = new BufferedReader(new FileReader(godl_freq_file));
            int num_sel_hap = 0;
            String g_line;
            while ((g_line = br_g.readLine()) != null) {
                g_line = g_line.replace("\n", "").replace("\r", "");
                final String[] tmp5 = g_line.split("\t");
                num_sel_hap = tmp5.length;
            }
            br_g.close();
            --num_sel_hap;
            final ArrayList<String> haps = new ArrayList<String>();
            final ArrayList<Double> haps_freq = new ArrayList<Double>();
            final ArrayList<String> tmp_haps = new ArrayList<String>();
            final ArrayList<Double> tmp_haps_freq = new ArrayList<Double>();
            final BufferedReader bufferedreader4 = new BufferedReader(new FileReader(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt"));
            final ArrayList<ArrayList<String>> geno0_2D = new ArrayList<ArrayList<String>>();
            String line3;
            while ((line3 = bufferedreader4.readLine()) != null) {
                line3 = line3.replace("\n", "").replace("\r", "");
                if (line3.startsWith("Hap_ID")) {
                    line3.split("\t");
                }
                if (line3.startsWith("Freq")) {
                    final String[] tmp6 = line3.split("\t");
                    for (int j = 1; j < tmp6.length; ++j) {
                        tmp_haps_freq.add(Double.parseDouble(tmp6[j]));
                    }
                }
                if (!line3.startsWith("Hap_ID") && !line3.startsWith("Freq")) {
                    final String[] tmp6 = line3.split("\t");
                    final ArrayList<String> tmp_arr3 = new ArrayList<String>();
                    for (int l = 1; l < tmp6.length; ++l) {
                        tmp_arr3.add(tmp6[l]);
                    }
                    geno0_2D.add(tmp_arr3);
                }
            }
            for (int m = 0; m < geno0_2D.get(0).size(); ++m) {
                String tmp_str = "";
                for (int l = 0; l < geno0_2D.size(); ++l) {
                    tmp_str = String.valueOf(tmp_str) + geno0_2D.get(l).get(m);
                }
                tmp_haps.add(tmp_str);
            }
            bufferedreader4.close();
            for (int i2 = 0; i2 < tmp_haps_freq.size(); ++i2) {
                for (int j2 = i2; j2 < tmp_haps_freq.size(); ++j2) {
                    if (tmp_haps_freq.get(i2) < tmp_haps_freq.get(j2)) {
                        final double tmp_freq = tmp_haps_freq.get(i2);
                        tmp_haps_freq.set(i2, tmp_haps_freq.get(j2));
                        tmp_haps_freq.set(j2, tmp_freq);
                        final String ss = tmp_haps.get(i2);
                        tmp_haps.set(i2, tmp_haps.get(j2));
                        tmp_haps.set(j2, ss);
                    }
                }
            }
            double total_freq = 1.0E-8;
            for (int l = 0; l < num_sel_hap; ++l) {
                if (l < tmp_haps_freq.size()) {
                    total_freq += tmp_haps_freq.get(l);
                    haps.add(tmp_haps.get(l));
                    haps_freq.add(tmp_haps_freq.get(l));
                }
            }
            for (int l = 0; l < haps_freq.size(); ++l) {
                haps_freq.set(l, haps_freq.get(l) / total_freq);
            }
            final BufferedWriter bw_gold = new BufferedWriter(new FileWriter(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt.snp10"));
            bw_gold.write("Hap_ID");
            for (int h = 0; h < haps.size(); ++h) {
                bw_gold.write("\th" + Integer.toString(h));
            }
            bw_gold.write("\nFreq");
            for (int h = 0; h < haps.size(); ++h) {
                bw_gold.write("\t" + haps_freq.get(h));
            }
            bw_gold.write("\n");
            for (int l2 = 0; l2 < haps.get(0).length(); ++l2) {
                bw_gold.write(index_var_prefix_dict2.get(l2));
                for (int h2 = 0; h2 < haps.size(); ++h2) {
                    bw_gold.write("\t" + haps.get(h2).substring(l2, l2 + 1));
                }
                bw_gold.write("\n");
            }
            bw_gold.close();
            haps.clear();
            haps_freq.clear();
            tmp_haps.clear();
            tmp_haps_freq.clear();
            final BufferedWriter bw_poolhapx_hap = new BufferedWriter(new FileWriter(String.valueOf(gp.out_dir) + "/poolhapx.haps"));
            if (index_var_prefix_dict2.size() > 50) {
                for (int pool_index15 = 0; pool_index15 < Entrance.num_pools; ++pool_index15) {
                    final ArrayList<ArrayList<String>> geno_2D = new ArrayList<ArrayList<String>>();
                    final String freq_out_file = String.valueOf(gp.out_dir) + "/" + Entrance.names_array[pool_index15] + "/final_freq_haps.txt";
                    System.out.println(freq_out_file);
                    final BufferedReader bufferedreader5 = new BufferedReader(new FileReader(freq_out_file));
                    boolean break_but = false;
                    while ((line3 = bufferedreader5.readLine()) != null) {
                        line3 = line3.replace("\n", "").replace("\r", "");
                        if (line3.startsWith("Freq")) {
                            break_but = false;
                            final String[] tmp7 = line3.split("\t");
                            for (int i3 = 1; i3 < tmp7.length; ++i3) {
                                if (Double.parseDouble(tmp7[i3]) > 1.0) {
                                    break_but = true;
                                }
                            }
                            if (break_but) {
                                bufferedreader5.close();
                                break;
                            }
                            for (int i3 = 1; i3 < tmp7.length; ++i3) {
                                tmp_haps_freq.add(Double.parseDouble(tmp7[i3]));
                            }
                        }
                        if (!line3.startsWith("Hap_ID") && !line3.startsWith("Freq")) {
                            final String[] tmp7 = line3.split("\t");
                            final ArrayList<String> tmp_arr4 = new ArrayList<String>();
                            for (int i4 = 1; i4 < tmp7.length; ++i4) {
                                tmp_arr4.add(tmp7[i4]);
                            }
                            geno_2D.add(tmp_arr4);
                        }
                    }
                    if (!break_but) {
                        for (int j3 = 0; j3 < geno_2D.get(0).size(); ++j3) {
                            String tmp_str2 = "";
                            for (int i4 = 0; i4 < geno_2D.size(); ++i4) {
                                tmp_str2 = String.valueOf(tmp_str2) + geno_2D.get(i4).get(j3);
                            }
                            tmp_haps.add(tmp_str2);
                        }
                    }
                    bufferedreader5.close();
                    for (int i5 = 0; i5 < tmp_haps_freq.size(); ++i5) {
                        for (int j4 = i5; j4 < tmp_haps_freq.size(); ++j4) {
                            if (tmp_haps_freq.get(i5) < tmp_haps_freq.get(j4)) {
                                final double tmp_freq2 = tmp_haps_freq.get(i5);
                                tmp_haps_freq.set(i5, tmp_haps_freq.get(j4));
                                tmp_haps_freq.set(j4, tmp_freq2);
                                final String ss2 = tmp_haps.get(i5);
                                tmp_haps.set(i5, tmp_haps.get(j4));
                                tmp_haps.set(j4, ss2);
                            }
                        }
                    }
                    total_freq = 1.0E-8;
                    for (int i5 = 0; i5 < default_sel_haps; ++i5) {
                        if (i5 < tmp_haps_freq.size()) {
                            total_freq += tmp_haps_freq.get(i5);
                            haps.add(tmp_haps.get(i5));
                            haps_freq.add(tmp_haps_freq.get(i5));
                        }
                    }
                    for (int i5 = 0; i5 < haps_freq.size(); ++i5) {
                        haps_freq.set(i5, haps_freq.get(i5) / total_freq);
                    }
                }
            }
            else {
                String freq_out_file2 = String.valueOf(gp.inter_dir) + "/aem/" + gp.project_name + "_level_3_region_0.inter_freq_haps.txt";
                if (index_var_prefix_dict2.size() > 20) {
                    freq_out_file2 = String.valueOf(gp.inter_dir) + "/aem/" + gp.project_name + "_level_5_region_0.inter_freq_haps.txt";
                }
                final BufferedReader bufferedreader6 = new BufferedReader(new FileReader(freq_out_file2));
                final ArrayList<ArrayList<String>> geno_2D2 = new ArrayList<ArrayList<String>>();
                while ((line3 = bufferedreader6.readLine()) != null) {
                    line3 = line3.replace("\n", "").replace("\r", "");
                    if (line3.startsWith("Freq")) {
                        final String[] tmp8 = line3.split("\t");
                        for (int i6 = 1; i6 < tmp8.length; ++i6) {
                            tmp_haps_freq.add(Double.parseDouble(tmp8[i6]));
                        }
                    }
                    if (!line3.startsWith("Hap_ID") && !line3.startsWith("Freq")) {
                        final String[] tmp8 = line3.split("\t");
                        final ArrayList<String> tmp_arr5 = new ArrayList<String>();
                        for (int i5 = 1; i5 < tmp8.length; ++i5) {
                            tmp_arr5.add(tmp8[i5]);
                        }
                        geno_2D2.add(tmp_arr5);
                    }
                }
                for (int j5 = 0; j5 < geno_2D2.get(0).size(); ++j5) {
                    String tmp_str3 = "";
                    for (int i5 = 0; i5 < geno_2D2.size(); ++i5) {
                        tmp_str3 = String.valueOf(tmp_str3) + geno_2D2.get(i5).get(j5);
                    }
                    tmp_haps.add(tmp_str3);
                }
                bufferedreader6.close();
                for (int i7 = 0; i7 < tmp_haps_freq.size(); ++i7) {
                    for (int j6 = i7; j6 < tmp_haps_freq.size(); ++j6) {
                        if (tmp_haps_freq.get(i7) < tmp_haps_freq.get(j6)) {
                            final double tmp_freq3 = tmp_haps_freq.get(i7);
                            tmp_haps_freq.set(i7, tmp_haps_freq.get(j6));
                            tmp_haps_freq.set(j6, tmp_freq3);
                            final String ss3 = tmp_haps.get(i7);
                            tmp_haps.set(i7, tmp_haps.get(j6));
                            tmp_haps.set(j6, ss3);
                        }
                    }
                }
                total_freq = 1.0E-8;
                for (int i7 = 0; i7 < default_sel_haps; ++i7) {
                    if (i7 < tmp_haps_freq.size()) {
                        total_freq += tmp_haps_freq.get(i7);
                        haps.add(tmp_haps.get(i7));
                        haps_freq.add(tmp_haps_freq.get(i7));
                    }
                }
                for (int i7 = 0; i7 < haps_freq.size(); ++i7) {
                    haps_freq.set(i7, haps_freq.get(i7) / total_freq);
                }
            }
            bw_poolhapx_hap.write("Hap_ID");
            for (int h2 = 0; h2 < haps.size(); ++h2) {
                bw_poolhapx_hap.write("\th" + Integer.toString(h2));
            }
            bw_poolhapx_hap.write("\nFreq");
            for (int h2 = 0; h2 < haps.size(); ++h2) {
                bw_poolhapx_hap.write("\t" + haps_freq.get(h2));
            }
            bw_poolhapx_hap.write("\n");
            for (int l3 = 0; l3 < haps.get(0).length(); ++l3) {
                bw_poolhapx_hap.write(index_var_prefix_dict2.get(l3));
                for (int h3 = 0; h3 < haps.size(); ++h3) {
                    bw_poolhapx_hap.write("\t" + haps.get(h3).substring(l3, l3 + 1));
                }
                bw_poolhapx_hap.write("\n");
            }
            bw_poolhapx_hap.close();
            eva.AEM_MCCEvaluate(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt.snp10", String.valueOf(gp.out_dir) + "/poolhapx.haps", freq_cutoff, String.valueOf(gp.out_dir) + "/poolhapx.MCC.result");
            eva.AEM_JSDEvaluate(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt.snp10", String.valueOf(gp.out_dir) + "/poolhapx.haps", 0.0, String.valueOf(gp.out_dir) + "/poolhapx.JSD.result");
            haps.clear();
            haps_freq.clear();
            tmp_haps.clear();
            tmp_haps_freq.clear();
            boolean do_aem = true;
            final File file4 = new File(String.valueOf(gp.out_dir) + "/aem.txt");
            if (!file4.exists()) {
                System.out.println(String.valueOf(gp.out_dir) + "/aem.txt");
                do_aem = false;
            }
            if (do_aem) {
                final BufferedReader bufferedreader7 = new BufferedReader(new FileReader(String.valueOf(gp.out_dir) + "/aem.txt"));
                while ((line3 = bufferedreader7.readLine()) != null) {
                    line3 = line3.replace("\n", "").replace("\r", "");
                    final String[] tmp8 = line3.split(" ");
                    String tmp_str3 = "";
                    if (!tmp8[tmp8.length - 1].equals("NA")) {
                        final double freq = Double.parseDouble(tmp8[tmp8.length - 1]);
                        if (freq <= 1.0E-9) {
                            continue;
                        }
                        for (int i4 = 0; i4 < tmp8.length - 1; ++i4) {
                            tmp_str3 = String.valueOf(tmp_str3) + tmp8[i4];
                        }
                        tmp_haps.add(tmp_str3);
                        tmp_haps_freq.add(freq);
                    }
                    else {
                        do_aem = false;
                    }
                }
                bufferedreader7.close();
            }
            if (do_aem) {
                for (int i8 = 0; i8 < tmp_haps_freq.size(); ++i8) {
                    for (int j5 = i8; j5 < tmp_haps_freq.size(); ++j5) {
                        if (tmp_haps_freq.get(i8) < tmp_haps_freq.get(j5)) {
                            final double tmp_freq4 = tmp_haps_freq.get(i8);
                            tmp_haps_freq.set(i8, tmp_haps_freq.get(j5));
                            tmp_haps_freq.set(j5, tmp_freq4);
                            final String ss4 = tmp_haps.get(i8);
                            tmp_haps.set(i8, tmp_haps.get(j5));
                            tmp_haps.set(j5, ss4);
                        }
                    }
                }
                total_freq = 1.0E-8;
                for (int i8 = 0; i8 < default_sel_haps; ++i8) {
                    if (i8 < tmp_haps_freq.size()) {
                        total_freq += tmp_haps_freq.get(i8);
                        haps.add(tmp_haps.get(i8));
                        haps_freq.add(tmp_haps_freq.get(i8));
                    }
                }
                for (int i8 = 0; i8 < haps_freq.size(); ++i8) {
                    haps_freq.set(i8, haps_freq.get(i8) / total_freq);
                }
                final BufferedWriter bw_aem_hap = new BufferedWriter(new FileWriter(String.valueOf(gp.out_dir) + "/aem.haps"));
                bw_aem_hap.write("Hap_ID");
                for (int h4 = 0; h4 < haps.size(); ++h4) {
                    bw_aem_hap.write("\th" + Integer.toString(h4));
                }
                bw_aem_hap.write("\nFreq");
                for (int h4 = 0; h4 < haps.size(); ++h4) {
                    bw_aem_hap.write("\t" + haps_freq.get(h4));
                }
                bw_aem_hap.write("\n");
                for (int l4 = 0; l4 < haps.get(0).length(); ++l4) {
                    bw_aem_hap.write(index_var_prefix_dict2.get(l4));
                    for (int h5 = 0; h5 < haps.size(); ++h5) {
                        bw_aem_hap.write("\t" + haps.get(h5).substring(l4, l4 + 1));
                    }
                    bw_aem_hap.write("\n");
                }
                bw_aem_hap.close();
                eva.AEM_MCCEvaluate(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt.snp10", String.valueOf(gp.out_dir) + "/aem.haps", freq_cutoff / 2.0, String.valueOf(gp.out_dir) + "/aem.MCC.result");
                eva.AEM_JSDEvaluate(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt.snp10", String.valueOf(gp.out_dir) + "/aem.haps", 0.005, String.valueOf(gp.out_dir) + "/aem.JSD.result");
            }
            haps.clear();
            haps_freq.clear();
            tmp_haps.clear();
            tmp_haps_freq.clear();
            final BufferedReader bufferedreader8 = new BufferedReader(new FileReader(String.valueOf(gp.out_dir) + "/hippo/results.out"));
            while ((line3 = bufferedreader8.readLine()) != null) {
                line3 = line3.replace("\n", "").replace("\r", "");
                final String[] tmp8 = line3.split(" ");
                final double freq2 = Double.parseDouble(tmp8[1]);
                if (freq2 > 1.0E-12) {
                    final String tmp_str2 = tmp8[0];
                    tmp_haps.add(tmp_str2);
                    tmp_haps_freq.add(freq2);
                }
            }
            for (int i7 = 0; i7 < tmp_haps_freq.size(); ++i7) {
                for (int j6 = i7; j6 < tmp_haps_freq.size(); ++j6) {
                    if (tmp_haps_freq.get(i7) < tmp_haps_freq.get(j6)) {
                        final double tmp_freq3 = tmp_haps_freq.get(i7);
                        tmp_haps_freq.set(i7, tmp_haps_freq.get(j6));
                        tmp_haps_freq.set(j6, tmp_freq3);
                        final String ss3 = tmp_haps.get(i7);
                        tmp_haps.set(i7, tmp_haps.get(j6));
                        tmp_haps.set(j6, ss3);
                    }
                }
            }
            for (int i7 = 0; i7 < default_sel_haps; ++i7) {
                if (i7 < tmp_haps_freq.size()) {
                    haps.add(tmp_haps.get(i7));
                    haps_freq.add(tmp_haps_freq.get(i7));
                }
            }
            bufferedreader8.close();
            final BufferedWriter bw_hippo_hap = new BufferedWriter(new FileWriter(String.valueOf(gp.out_dir) + "/hippo.haps"));
            bw_hippo_hap.write("Hap_ID");
            for (int h5 = 0; h5 < haps.size(); ++h5) {
                bw_hippo_hap.write("\th" + Integer.toString(h5));
            }
            bw_hippo_hap.write("\nFreq");
            for (int h5 = 0; h5 < haps.size(); ++h5) {
                bw_hippo_hap.write("\t" + haps_freq.get(h5));
            }
            bw_hippo_hap.write("\n");
            for (int l5 = 0; l5 < haps.get(0).length(); ++l5) {
                bw_hippo_hap.write(index_var_prefix_dict2.get(l5));
                for (int h6 = 0; h6 < haps.size(); ++h6) {
                    bw_hippo_hap.write("\t" + haps.get(h6).substring(l5, l5 + 1));
                }
                bw_hippo_hap.write("\n");
            }
            bw_hippo_hap.close();
            eva.AEM_MCCEvaluate(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt.snp10", String.valueOf(gp.out_dir) + "/hippo.haps", 0.0, String.valueOf(gp.out_dir) + "/hippo.MCC.result");
            eva.AEM_JSDEvaluate(String.valueOf(gp.gold_dir) + "/" + gp.project_name + "_haps.inter_freq_vars.txt", String.valueOf(gp.out_dir) + "/hippo.haps", 0.0, String.valueOf(gp.out_dir) + "/hippo.JSD.result");
        }
    }
}