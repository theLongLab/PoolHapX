package Viral_Reconstructions_Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class Transfer_output_to_fasta {
	/**
	 *  This class is used to transfer the output into fasta file in order for
	 *  using the clustalo to do the mapping
	 *  
	 *  The final output format looks like this:
	 *  >Reference
	 *  TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACA
	 *  >Hap0_pool0_0.4
	 *  TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACA
	 *  >Hap1_pool0_0.6
	 *  TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACA
	 *  >Hap0_pool1_0.5
	 *  TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACA
	 *  >Hap1_pool1_0.5
	 *  TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACA
	**/
	
	  String ref_seq_file_path;
	  String output_dir;
	  String project_name;
	  String output_suffix; 
	  String fasta_folder;
	  int num_pools;
	  
	public Transfer_output_to_fasta(String parameter_file) throws IOException {
	    InputStream is = new FileInputStream(parameter_file);
		Properties prop = new Properties();
		prop.load(is);
		this.output_dir = prop.getProperty("Output_Dir");
		this.ref_seq_file_path = prop.getProperty("Ref_Seq_Path");
		this.output_suffix = prop.getProperty("Output_Suffix");		
		this.project_name=prop.getProperty("Proj_Name");
		this.num_pools= Integer.parseInt(prop.getProperty("Num_Pools"));
		this.fasta_folder=this.output_dir+"fasta/";
		 is.close();
		}
		
//	public static void output_to_fasta(String ref_seq_file_path, 
//			String output_dir, String project_name, String output_suffix, 
//			int num_pools) throws IOException, InterruptedException {
		// First read the reference_sequence file
	public void output_to_fasta() throws IOException, InterruptedException{
		BufferedReader br_ref_file = new BufferedReader(new FileReader(
				ref_seq_file_path));
		BufferedWriter bw_fasta_file = new BufferedWriter(new FileWriter(
				fasta_folder + project_name + ".fasta"));
		    bw_fasta_file.write(">Reference"+"\n");
			String ref_line=br_ref_file.readLine();
			ref_line=br_ref_file.readLine();
			while(ref_line!=null) {
				String [] each_position=ref_line.split("");
				for(int i =0 ;i < each_position.length; i++) {
					bw_fasta_file.write(each_position[i]);
				}
				ref_line=br_ref_file.readLine();
			}
			br_ref_file.close();
        
		for (int p=0; p < num_pools; p++) {
			String sample_names = project_name + "_p" + p; // 0_1_p0
			String output_file = sample_names + "." + output_suffix; // 0_1_p0.fasta
			int num_haps_per_pool = 0;
			BufferedReader br = new BufferedReader(new FileReader(output_dir 
					+ output_file));
			String currline = br.readLine(); // read the first line
			while(currline!=null) {
				
				String[] each_position = currline.split("");
				if(each_position[0].equals(">")) {
					bw_fasta_file.write("\n");
					String[] strain_fre_line=currline.split("_");
					bw_fasta_file.write(">"+ "Hap" + num_haps_per_pool
							+ "_pool"+ p + "_" 
					+ strain_fre_line[strain_fre_line.length-1] + "\n" );
					num_haps_per_pool++;
				}else {
					for(int i =0; i < each_position.length; i++) {
						bw_fasta_file.write(each_position[i]);
					}
					
				}
				currline = br.readLine();
			}
			br.close();
		}	
		bw_fasta_file.close();
	}
	
	public static void main(String[] args) throws IOException, InterruptedException {
		String parameter = args[0];//"D:\\PhD-Studying\\Informatics\\Project\\HIV project\\Viral_reconstruction\\Other_tools_results\\CliqueSNV\\fasta\\FS3.properties";
		Transfer_output_to_fasta gf = new Transfer_output_to_fasta(parameter);
		gf.output_to_fasta();
	}
	
	
}



