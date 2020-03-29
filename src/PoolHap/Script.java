package PoolHap;

import java.io.Writer;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.Reader;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Properties;
import java.io.FileInputStream;
import java.util.ArrayList;

public class Script
{
    public String main_dir;
    public String project_name;
    public String java;
    public String phx_jar;
    public String bwa;
    public String samtools;
    public String gatk;
    public String fastq_path;
    public String fastq_file;
    public String ref_path;
    public String Rscript_path;
    public String species;
    public ArrayList<String> sample_name_list;
    public String seq_tech; 
    public String longranger; 
    public String longranger_ref_folder;
    
    public Script(final String config_path) throws IOException {
        this.sample_name_list = new ArrayList<String>();
        final InputStream is = new FileInputStream(config_path);
        final Properties prop = new Properties();
        prop.load(is);
        this.main_dir = prop.getProperty("Main_Dir");
        this.project_name = prop.getProperty("Project_Name");
        this.java = prop.getProperty("Java");
        this.phx_jar = prop.getProperty("PHX_JAR");
        this.bwa = prop.getProperty("bwa");
        this.longranger= prop.getProperty("longranger");
        this.samtools = prop.getProperty("samtools");
        this.gatk = prop.getProperty("gatk");
        this.fastq_path = prop.getProperty("Fastq_Path");
        this.ref_path = prop.getProperty("Ref_Path");
        this.longranger_ref_folder = prop.getProperty("Longranger_Ref_Folder");
        this.Rscript_path = prop.getProperty("Rscript_Path");
        this.fastq_file = prop.getProperty("Fastq_File");
        this.seq_tech= prop.getProperty("Sequencing_Technology");
        is.close();
//        System.out.println(this.fastq_file);
//        System.out.println(this.fastq_file);
        final BufferedReader br = new BufferedReader(new FileReader(this.fastq_file));
        for (String currline = br.readLine(); currline != null; currline = br.readLine()) {
        	if (this.seq_tech.equals("paired-end_reads")) {
	            final String[] tmpcurrpos = currline.split("\t");
	            final String curr_sample_name = tmpcurrpos[0].split(".read1.fastq")[0];
	            this.sample_name_list.add(curr_sample_name);
        	}
        	if (this.seq_tech.equals("10x_linked-reads")) {
        		this.sample_name_list.add(currline);
        	}
        }
        br.close();
    }
    
    public void ComparisonSim() throws IOException {
        new File(String.valueOf(this.main_dir) + "/cmd/").mkdir();
        
        final String project_name = this.project_name;
        new File(String.valueOf(this.main_dir) + "/" + project_name).mkdir();
        new File(String.valueOf(this.main_dir) + "/" + project_name + "/input/").mkdir();
        
        final BufferedWriter bw = new BufferedWriter(new FileWriter(String.valueOf(this.main_dir) + "/cmd/" + project_name + ".cmd"));
        
        bw.write("java=" + this.java + "\n");
        bw.write("poolhapx=" + this.phx_jar + "\n");
        bw.write("bwa=" + this.bwa + "\n");
        bw.write("samtools=" + this.samtools + "\n");
        bw.write("gatk=" + this.gatk + "\n");
        bw.write("ref=" + this.ref_path + "\n");
        bw.write("longranger=" + this.longranger + "\n"); 
        bw.write("longranger_ref_folder=" + this.longranger_ref_folder + "\n"); 
        bw.write("mkdir " + this.main_dir + "/" + project_name + "/input/bam\n");
        bw.write("mkdir " + this.main_dir + "/" + project_name + "/input/vcf\n");
        bw.write("mkdir " + this.main_dir + "/" + project_name + "/input/sam\n");
        bw.write("mkdir " + this.main_dir + "/" + project_name + "/output\n");
        bw.write("mkdir " + this.main_dir + "/" + project_name + "/intermediate\n");
        for (int i = 0; i < this.sample_name_list.size(); ++i) {
            final String sample_name = this.sample_name_list.get(i);
            bw.write("prefix_fastq=" + this.fastq_path + "/" + sample_name + "\n");
            bw.write("prefix_bam=" + this.main_dir + "/" + project_name + "/input/bam/" + sample_name + "\n");
            bw.write("inbam=" + this.main_dir + "/" + project_name + "/input/bam/" + sample_name + ".rg.bam\n");
            bw.write("prefix_vcf=" + this.main_dir + "/" + project_name + "/input/vcf/" + sample_name + "\n");
            bw.write("outgvcf=$prefix_vcf\\.raw.g.vcf\n");
//            /export/qlong/PoolHapX/longranger-2.2.2/longranger  align --id=freq_4_pool_50_dep_250_p4  --reference=/export/qlong/chencao/Work/poolhapx/slim/ref/plasmodium/refdata-plasmodium  --fastqs=$prefix
            if (this.seq_tech.equals("10x_linked-reads")) {
            	bw.write("cd $prefix_fastq\n");
            	bw.write("$longranger align --id="+sample_name+" --reference=$longranger_ref_folder --fastqs=$prefix_fastq\n");            	
            	bw.write("$samtools sort -o  $prefix_bam\\.srt.bam  $prefix_fastq"+"\\/"+ sample_name + "/outs/possorted_bam.bam\n");
            	bw.write("$gatk  AddOrReplaceReadGroups -I  $prefix_bam\\.srt.bam -O $inbam -R $ref -ID " + sample_name + " -LB NPD -PL Illumina -PU NPD -SM " + sample_name + "\n");
            }
            
            if (this.seq_tech.equals("paired-end_reads")) {
            	bw.write("$bwa mem $ref $prefix_fastq\\.read1.fastq $prefix_fastq\\.read2.fastq | $samtools view -Shub - > $prefix_bam\\.bam\n");
            	bw.write("$samtools sort -o  $prefix_bam\\.srt.bam  $prefix_bam\\.bam\n");
            	bw.write("$gatk  AddOrReplaceReadGroups -I  $prefix_bam\\.srt.bam -O $inbam -R $ref -ID " + sample_name + " -LB NPD -PL Illumina -PU NPD -SM " + sample_name + "\n");
            }

            
            bw.write("$samtools index $inbam\n");
            bw.write("$gatk  HaplotypeCaller -R  $ref -I $inbam -ERC  GVCF -ploidy 8 --heterozygosity 0.01  --max-alternate-alleles 1 -O $outgvcf\n");
        }
        String tmp = "$gatk CombineGVCFs -R $ref ";
        for (int j = 0; j < this.sample_name_list.size(); ++j) {
            final String sample_name2 = this.sample_name_list.get(j);
            tmp = String.valueOf(tmp) + " -V " + this.main_dir + "/" + project_name + "/input/vcf/" + sample_name2 + ".raw.g.vcf";
        }
        bw.write("prefix_project_vcf=" + this.main_dir + "/" + project_name + "/input/vcf/" + project_name + "\n");
        tmp = String.valueOf(tmp) + " -O $prefix_project_vcf\\.g.vcf\n";
        bw.write(tmp);
        bw.write("$gatk GenotypeGVCFs -R $ref -V  $prefix_project_vcf\\.g.vcf -ploidy 8 -O $prefix_project_vcf\\.raw.vcf\n");
        bw.write("prefix_project=" + this.main_dir + "/" + project_name + "/input/" + project_name + "\n");
        bw.write("$gatk  SelectVariants -R $ref -V  $prefix_project_vcf\\.raw.vcf -O $prefix_project\\.vcf\n");
        for (int j = 0; j < this.sample_name_list.size(); ++j) {
            final String sample_name2 = this.sample_name_list.get(j);
            bw.write("prefix_sam=" + this.main_dir + "/" + project_name + "/input/sam/" + sample_name2 + "\n");
            bw.write("prefix_bam=" + this.main_dir + "/" + project_name + "/input/bam/" + sample_name2 + "\n");
            bw.write("$samtools view -ho $prefix_sam\\.sam $prefix_bam\\.srt.bam\n");
        }
        final BufferedWriter bw_properties = new BufferedWriter(new FileWriter(String.valueOf(this.main_dir) + "/" + project_name + "/input/PHX.properties"));
        bw_properties.write("Proj_Name = " + project_name + "\n");
        bw_properties.write("Input_Dir = " + this.main_dir + "/" + project_name + "/input/\n");
        bw_properties.write("Intermediate_Dir = " + this.main_dir + "/" + project_name + "/intermediate/\n");
        bw_properties.write("Output_Dir = " + this.main_dir + "/" + project_name + "/output/\n");
        bw_properties.write("Num_Pos_Window = 20\n");
        bw_properties.write("Num_Gap_Window = 2\n");
        bw_properties.write("In-pool_Gap_Support_Min = 1\n");
        bw_properties.write("All-pool_Gap_Support_Min = 1\n");
        bw_properties.write("Level_1_Region_Size_Min = 10\n");
        bw_properties.write("Level_1_Region_Size_Max = 12\n");
        bw_properties.write("Level_1T_Region_Size_Min = 10\n");
        bw_properties.write("Level_1T_Region_Size_Max = 12\n");
        bw_properties.write("Est_Ind_PerPool = 1000000\n");
        bw_properties.write("Level_1_2_Region_Mismatch_Tolerance = 1\n");
        bw_properties.write("Level_2_3_Region_Mismatch_Tolerance = 2\n");
        bw_properties.write("Level_3_4_Region_Mismatch_Tolerance = 5\n");
        bw_properties.write("AEM_Maximum_Level = 4\n");
        bw_properties.write("BFS_Mismatch_Tolerance = 6\n");
        bw_properties.write("AEM_Iterations_Max = 200\n");
        bw_properties.write("AEM_Convergence_Cutoff = 0.00001\n");
        bw_properties.write("AEM_Zero_Cutoff = 0.00001\n");
        bw_properties.write("AEM_Regional_Cross_Pool_Freq_Cutoff = 0.01\n");
        bw_properties.write("AEM_Regional_HapSetSize_Max = 50\n");
        bw_properties.write("AEM_Regional_HapSetSize_Min = 3\n");
        bw_properties.write("IF_0_0 = 0.1\n");
        bw_properties.write("IF_Denominator_0= 10.0\n");
        bw_properties.write("Rscript_path =" + this.Rscript_path + "\n");
        bw_properties.write("Regression_Distance_Max_Weight = 2.0\n");
        bw_properties.write("Regression_Coverage_Weight = 2.0\n");
        bw_properties.write("Regression_One_Vector_Weight = 5.0\n");
        bw_properties.write("Regression_Hap_MAF_Weight = 2.0\n");
        bw_properties.write("Regression_Hap_LD_Weight = 1.0\n");
        bw_properties.write("Regression_Mismatch_Tolerance =7\n");
        bw_properties.write("Maximum_Selected_HapSetSize = 25\n");
        bw_properties.write("Regression_Gamma_Min = 0.0001\n");
        bw_properties.write("Regression_Gamma_Max = 0.1\n");
        bw_properties.write("Regression_n_Gamma = 10\n");
        bw_properties.write("Regression_Maximum_Regions = 3\n");
        
        if (this.seq_tech.equals("paired-end_reads")) {
        	bw_properties.write("Sequencing_Technology=paired-end_reads\n");
        }
        if (this.seq_tech.equals("10x_linked-reads")) {
        	bw_properties.write("Sequencing_Technology=10x_linked-reads\n");
        }
        
        bw_properties.write("Number_Threads = 3\n");
        bw_properties.close();
        bw.write("properties=" + this.main_dir + "/" + project_name + "/input/PHX.properties\n");
        bw.write("poolhapx=" + this.phx_jar + "\n");
        bw.write("$java -jar $poolhapx format $properties\n");
        bw.write("$java -jar $poolhapx gc $properties\n");
        bw.write("$java -jar $poolhapx aem $properties\n");
        bw.write("$java -jar $poolhapx l0l1 $properties\n");
        bw.close();
        
        Process ps = Runtime.getRuntime().exec("chmod 777 "+String.valueOf(this.main_dir) + "/cmd/" + 
        		project_name + ".cmd");
    }
}