package MiscFunctions;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

public class BAMFormatterAlt {

    public static void main(String[] args) throws IOException {

        String outdir = args[0]; // inter_dir/
        String prefix = args[1]; // combination_simulation
        Integer numPts = Integer.parseInt(args[2]);	// number of pools (patients) // TODO: rename?

        // TODO: [Question]:: what is this? leftover? ML 20190702
        Boolean simmode = false; // Boolean.parseBoolean(args[3]);

        HashSet<Integer> trueVarPos = new HashSet<Integer>();

        // TODO: LEFTOVER(?) ML 20190702
        if (simmode) {
            trueVarPos = simvarsReader(outdir + "/simvars.intra_freq.txt", numPts);
        }

        PrintWriter VEFList = new PrintWriter(outdir + prefix + ".vef.list");

        // Encode variants from given VCF file.
        HashMap<Integer, HashMap<String, VarObj>> variantEncoder = VariantMapping(
            // TODO: [I/O]:: remove filename hardcoding. ML 20190702
            outdir + prefix + ".vcf", // args[1] needs to be interdir + prefix + .all.vcf.
            true,
            numPts,
            simmode,
            trueVarPos);

        int numVars = variantEncoder.size();
        double[][] poolFreqs = new double[numVars][numPts];
        for (int p = 0; p < numPts; p++) { // set pool frequency file
            // TODO: [I/O]:: remove filename hardcoding. ML 20190702
            double[] pool_freq = VEFMaker(outdir + prefix + "_p" + p, variantEncoder);
            VEFList.println(outdir + prefix + "_p" + p + ".vef");

            for (int v = 0; v < numVars; v++) {
                poolFreqs[v][p] = pool_freq[v];
            }
        }
        VEFList.close();

        BufferedWriter br = new BufferedWriter(new FileWriter(outdir
            + prefix
            + "_vars.intra_freq.txt"));

        // Intra-pool variant frequency file header.
        br.write("Pool_ID\t");
        for (int p = 0; p < numPts; p++) {
            br.write(p + "\t");
        }

        ArrayList<Integer> posTracker = new ArrayList<Integer>();
        posTracker.addAll(variantEncoder.keySet()); // add variant positions
        Collections.sort(posTracker); // sort variant positions
        for (int v = 0; v < numVars; v++) { // for each variant position...
            // TODO: [old] This only allows for biallelic simple loci (single alternate allele) for
            // now.
            for (int p = 0; p < numPts; p++) { // for each pool, write the variant frequency
                br.write("\t" + poolFreqs[v][p]);
            }

            // TODO: LEFTOVER ML 20190702
            // br.write("\n0;" + posTracker.get(v) + ";" + posTracker.get(v) + ";0:1");

        }
        br.close();

    }

    // TODO: LEFTOVER(?) ML 20190702
    private static HashSet<Integer> simvarsReader(
        String simvarsFile,
        int numPts) throws IOException {

        BufferedReader SVReader = new BufferedReader(new FileReader(simvarsFile));
        String currLine = SVReader.readLine();	// Skip the first line.
        currLine = SVReader.readLine();
        HashSet<Integer> trueVarPos = new HashSet<Integer>();
        svRead: while (currLine != null) {
            String[] fullLine = currLine.split("\t");
            // 0;211;211;0:1		0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
            for (int p = 2; p <= numPts + 1; p++) {
                if (Double.parseDouble(fullLine[p]) == 0)  {
                    currLine = SVReader.readLine();
                    continue svRead;
                }
            }
            int truePos = Integer.parseInt(fullLine[0].split(";")[1]);
            trueVarPos.add(truePos);
            currLine = SVReader.readLine();
        }
        SVReader.close();
        return trueVarPos;
    }

    /**
     *  Encode variants within VCF file into dictionary of position to alternate allele.
     *
     *  @param VCFFile
     *  @param biallelic
     *  @param numPts
     *  @param simmode
     *  @param trueVarPos
     *  @return
     *  @throws IOException
     */
    private static HashMap<Integer, HashMap<String, VarObj>> VariantMapping(
        String VCFFile,
        boolean biallelic,
        int numPts,
        Boolean simmode,
        HashSet<Integer> trueVarPos) throws IOException {

        BufferedReader VCFReader = new BufferedReader(new FileReader(VCFFile));
        String currLine = VCFReader.readLine();

        // TODO: LEFTOVER ML 20190702
        // System.out.println(currLine);

        // Skip VCF metadata lines.
        while (!currLine.contains("#CHROM")) {
            currLine = VCFReader.readLine();
        }

        // Match up the order of the pool IDs in the VCF to their actual row number.
        String[] header_line = currLine.split("\t"); // store header
        int[] actual_row = new int[numPts];

        // Read pool names. // TODO: [ReconEP] make it so names don't have to be integers.
        // ML 20190702
        for (int i = 9; i < header_line.length; i++) {
            actual_row[i - 9] = Integer.parseInt(header_line[i]);
        }

        currLine = VCFReader.readLine(); // first data line

        // TODO: LEFTOVER ML 20190702
        // System.out.println(currLine);

        HashMap<Integer,HashMap<String,VarObj>> variantEncoder =
            new HashMap<Integer, HashMap<String, VarObj>>();

        // TODO: LEFTOVER ML 20190702
        // ArrayList<double[]> poolFreqs = new ArrayList<double[]>();
        // ArrayList<Integer> posTracker = new ArrayList<Integer>();
        // int v = 0;

        // Read through rest of VCF file.
        while (currLine != null) {
            String[] var_info = currLine.split("\t");
            int pos = Integer.parseInt(var_info[1]);

            // TODO: LEFTOVER(?) ML 20190702
            if (simmode) {
                if (!trueVarPos.contains(pos)) {
                    currLine = VCFReader.readLine();
                    continue;	// If this is not a true variant position, skip it.
                }
            }

            // TODO: LEFTOVER ML 20190702
            // posTracker.add(pos);
            // System.out.println(pos);

            variantEncoder.put(pos, new HashMap<String, VarObj>());

            // This Scanner runs through the String containing all of the alternate allele(s).
            // There is at least one.
            Scanner altAlleleScanner = new Scanner(var_info[4]);
            altAlleleScanner.useDelimiter(",");
            int variantCode = 1; // always going to be 1 because PHX discards the rest at the moment

            // In case of highly polymorphic (>=2 alternate alleles) positions.
            while (altAlleleScanner.hasNext()) {
                String variantNow = altAlleleScanner.next().replace("\"","");

                // This is the information that will be reported when the alternate allele at this
                // position in the read is accessed.
                VarObj altAlleleAtPos = new VarObj(pos, variantCode);

                // The information is mapped to the alternate allele.
                variantEncoder.get(pos).put(variantNow, altAlleleAtPos);

                // The code corresponding to the alternate allele is incremented so no two
                // alternate alleles in the same location have the same code.
                variantCode++;

                if (biallelic == true) {
                    break; 	// If PoolHap has not been adjusted for 3+ alleles.
                }

            }
            altAlleleScanner.close();

            // TODO: LEFTOVER ML 20190702
            // poolFreqs.add(new double[numPts]);
            // for (int p = 0; p < numPts; p++) {
            //     // This is the AD of the GT:AD:DP:GQ:PL of the p0 block of genotypes.
            //     String[] varCts = var_info[p + 9].split(":")[1].split(",");

            //     double ref = Double.parseDouble(varCts[0]);
            //     double alt = Double.parseDouble(varCts[1]);
            //     if (ref != 0 || alt != 0) {
            //         poolFreqs.get(v)[actual_row[p]] = alt / (ref + alt);

            //     // In case there are no reads at this position at all.
            //     } else {
            //         poolFreqs.get(v)[actual_row[p]] = 0;
            //     }
            //     System.out.println(poolFreqs.get(v)[p]);
            // }
            // v++;

            currLine = VCFReader.readLine();

        }
        VCFReader.close();

        // TODO: LEFTOVER ML 20190702
        // System.out.println("Done.");

        return variantEncoder;

    }


    /**
     *
     *
     *  @param BAMPrefix
     *  @param variantEncoder
     *  @return
     *  @throws FileNotFoundException
     */
    private static double[] VEFMaker(
        String BAMPrefix,
        HashMap<Integer, HashMap<String, VarObj>> variantEncoder) throws FileNotFoundException {

        Set<Integer> keyMATCH, keyINS, keyDEL, variantKS =
            variantEncoder.keySet(),
            tempSet,
            indelKS;

        SortedSet<Integer> variantSS = new TreeSet<Integer>();
        variantSS.addAll(variantKS);
        int numVars = variantEncoder.size();
        double[] pool_ct = new double[numVars];
        double[] total_ct = new double[numVars];
        File BAMPath = new File(BAMPrefix + ".srt.sam");
        Scanner BAMScanner = new Scanner(BAMPath);
        BAMScanner.useDelimiter("\t");
        String currLine = BAMScanner.nextLine();

        // TODO: LEFTOVER ML 20190702
        // System.out.println(currLine);

        // Skip headers.
        while (currLine.matches("@(.*)")) {
            String prevLine = currLine;
            currLine = BAMScanner.nextLine();
            if (prevLine.contains("PN:")) {
                break;
            }

            // TODO: LEFTOVER ML 20190702
            // System.out.println(currLine);

        }

        String QNAME, CIGAR, SEQ, tempPos = "", currAltAllele;
        Integer FLAG, startPOS = 1, currentPos, posToAdd, matchToIndel, skipSBases = 0;
        VarObj currVarObj;

        // hm = HashMap.
        HashMap<Integer, HashMap<Integer, String>> hmMATCH =
            new HashMap<Integer, HashMap<Integer, String>>();

        HashMap<Integer, HashMap<Integer, String>> hmINS =
            new HashMap<Integer, HashMap<Integer, String>>();

        HashMap<Integer, HashMap<Integer, String>> hmDEL =
            new HashMap<Integer, HashMap<Integer, String>>();

        HashMap<Integer, String> defaultHM, tempFinderHM;
        HashMap<String, VarObj> tempReporterHM;

        // TODO: LEFTOVER ML 20190702
        // ArrayList<double[]> poolFreqs = new ArrayList<double[]>();
        // HashMap<Integer,Integer> hmVarCount = new HashMap<Integer,Integer>();
        // for (Integer i : variantSS) {
        //     hmVarCount.put(i,0);
        // }
        // int readCount = 0;

        PrintWriter VEFFile = new PrintWriter(BAMPrefix + ".raw.vef");
        while (BAMScanner.hasNextLine()) {

            // TODO: LEFTOVER ML 20190702
            // readCount++;

            QNAME = BAMScanner.next().replace(":", "");

            // TODO: LEFTOVER ML 20190702
            // System.out.println(QNAME);

            FLAG = BAMScanner.nextInt();
            if (FLAG > 163) {
                // If the read is chimeric, or mapped to more than one location, or of poor
                // technical quality, skip it.
                BAMScanner.nextLine();
                continue;

            }
            BAMScanner.next(); // Skip the FLAG, RNAME.
            startPOS = BAMScanner.nextInt();
            currentPos = startPOS;
            BAMScanner.next(); // Skip the MAPQ.
            CIGAR = BAMScanner.next();

            if (CIGAR.charAt(0) == '*') {
                // * refers to the fact that no CIGAR information is available.
                BAMScanner.nextLine();

                // The above line skips the rest of the information in this read entirely.
                continue;
            }

            for (int c = 0; c < CIGAR.length(); c++) {
                Character tempChar = CIGAR.charAt(c);
                if (tempChar.compareTo('M') == 0) { // M = match
                    posToAdd = Integer.parseInt(tempPos);
                    for (int addPos = currentPos; addPos < currentPos + posToAdd; addPos++) {
                        defaultHM = new HashMap<Integer, String>();
                        defaultHM.put(1, "N");
                        hmMATCH.put(addPos, defaultHM);
                    }
                    currentPos += posToAdd;
                    tempPos = "";

                } else if (tempChar.compareTo('I') == 0) {
                    // Command to check for insertions:
                    // samtools view HIV_p1_sample_1.procd.bam | awk '{if ($6 ~ /I/) print $6;}' -
                    posToAdd = Integer.parseInt(tempPos);

                    // This is the first position of the insertion, which it is recognized by.
                    matchToIndel = currentPos - 1;

                    // Remove the first position of the insertion from the match HM.
                    hmMATCH.remove(matchToIndel);
                    defaultHM = new HashMap<Integer, String>();

                    // This accounts for the entire length of the insertion.
                    defaultHM.put(posToAdd + 1, "N");
                    hmINS.put(matchToIndel, defaultHM);

                    // Since the insertion starts at the previous currentPos - 1, the next match
                    // along starts at matchToIndel + 1 = previous currentPos.
                    currentPos = matchToIndel + 1;
                    tempPos = "";

                } else if (tempChar.compareTo('D') == 0) {
                     // This is the first position of the deletion, which it is recognized by.
                    posToAdd = Integer.parseInt(tempPos);
                    matchToIndel = currentPos - 1;

                    // Remove the first position of the deletion from the match HM.
                    hmMATCH.remove(matchToIndel);
                    defaultHM = new HashMap<Integer, String>();

                    // This accounts for the entire length of the deletion.
                    defaultHM.put(posToAdd + 1, "N");
                    hmDEL.put(matchToIndel, defaultHM);

                    // Since the deletion starts at the previous currentPos - 1, the next match
                    // along starts at currentPos + matchToIndel - 1.
                    currentPos += posToAdd;
                    tempPos = "";

                } else if (tempChar.compareTo('S') == 0) {
                    // These bases do appear in SEQ BUT they are not aligned to the reference
                    // genome. Therefore these bases do not contribute to the currentPOS.
                    skipSBases = Integer.parseInt(tempPos);
                    tempPos = "";

                } else if (tempChar.compareTo('H') == 0) {
                    // These bases DO NOT appear in SEQ and they are not aligned to the reference
                    // genome. Therefore these bases do not contribute to the currentPOS.
                    tempPos = "";

                } else {
                    tempPos += tempChar;
                }
            }

            keyMATCH = hmMATCH.keySet();
            keyINS = hmINS.keySet();
            keyDEL = hmDEL.keySet();
            for (int s = 0; s < 3; s++) {
                BAMScanner.next(); 	// Skip the RNEXT, PNEXT, TLEN.
            }
            SEQ = BAMScanner.next();
            int endPOS = startPOS + SEQ.length() - 1;
            currentPos = startPOS;

            // 4. Test the following to see if the SEQ-parsing code works properly, and the proper
            // bases are assigned to the correct Hash Maps.
            for (int s = skipSBases; s < SEQ.length(); s++) {
                // Skip all of the bases that were soft-clipped but still reported at the beginning
                // of SEQ. In the event there is no 'S' in the CIGAR string, this is 0 and we start
                // reporting from the end of SEQ.

                // TODO: LEFTOVER ML 20190702
                // System.out.println(s + "\t" + currentPos);

                if (keyMATCH.contains(currentPos)) {
                    hmMATCH.get(currentPos).put(1, SEQ.substring(s, s + 1));
                    currentPos++;

                } else if (keyINS.contains(currentPos)) {
                    tempSet = hmINS.get(currentPos).keySet();
                    for (Integer length : tempSet) {
                        hmINS.get(currentPos).put(length, SEQ.substring(s, s + length));
                    }
                    currentPos++;

                } else if (keyDEL.contains(currentPos)) {
                    tempSet = hmDEL.get(currentPos).keySet();
                    for (Integer length : tempSet) {
                        hmDEL.get(currentPos).put(length, SEQ.substring(s, s + 1));
                        currentPos += length;
                    }
                }
            }

            // TODO: LEFTOVER ML 20190702
            // for (int i : keyMATCH) {
            //     System.out.println(i + "\t" + hmMATCH.get(i));
            // }

            // Reset the 'start' position of SEQ-processing to 0 in case there aren't any S-es in
            // the next CIGAR string.
            skipSBases = 0;

            StringBuilder readInfo = new StringBuilder();

            // TODO: LEFTOVER ML 20190702
            // System.out.println(startPOS + "\t" + endPOS);

            int v = 0;
            // subSet method is (inclusive, exclusive)
            for (Integer VCFVarPos : variantSS.subSet(startPOS, endPOS + 1)) {
                // TODO: LEFTOVER ML 20190702
                // System.out.println(VCFVarPos);

                tempReporterHM = variantEncoder.get(VCFVarPos);
                if (keyMATCH.contains(VCFVarPos)) {
                    currAltAllele = hmMATCH.get(VCFVarPos).get(1);

                    // TODO: LEFTOVER ML 20190702
                    // System.out.println(currAltAllele);

                    if (tempReporterHM.containsKey(currAltAllele)) {
                        currVarObj = tempReporterHM.get(currAltAllele);

                        // TODO: LEFTOVER ML 20190702
                        // VEFFile.print(currVarObj.intM + "=" + currVarObj.varCode + ";");

                        readInfo.append(currVarObj.intM + "=" + currVarObj.varCode + ";");

                        // TODO: LEFTOVER ML 20190702
                        // hmVarCount.put(VCFVarPos, hmVarCount.get(VCFVarPos) + 1);

                        pool_ct[v]++;

                    } else {
                        // Fixed as of 10282018. Noted by CC that in some reads, variant positions
                        // in those reads were not being annotated. Realized that the 'else' clause
                        // was in the wrong place i.e.: alleles  =/= alternate (i.e: the reference)
                        // were not noted when they existed.
                        readInfo.append(VCFVarPos + "=0;");
                        // System.out.println(VCFVarPos + " done");
                    }
                    total_ct[v]++;

                    // TODO: LEFTOVER ML 20190702
                    // System.out.println("match");

                } else if (keyINS.contains(VCFVarPos)) {
                    tempFinderHM = hmINS.get(VCFVarPos);
                    indelKS = tempFinderHM.keySet();
                    for (Integer indel : indelKS) {
                        currAltAllele = tempFinderHM.get(indel);
                        if (tempReporterHM.containsKey(currAltAllele)) {
                            currVarObj = tempReporterHM.get(currAltAllele);

                            // TODO: LEFTOVER ML 20190702
                            // VEFFile.print(currVarObj.intM + "=" + currVarObj.varCode + ";");

                            readInfo.append(currVarObj.intM + "=" + currVarObj.varCode + ";");

                            // TODO: LEFTOVER ML 20190702
                            // hmVarCount.put(VCFVarPos, hmVarCount.get(VCFVarPos) + 1);

                            pool_ct[v]++;

                        } else {
                            readInfo.append(VCFVarPos + "=0;");

                            // TODO: LEFTOVER ML 20190702
                            // System.out.println(VCFVarPos + " done");

                        }
                    }
                    total_ct[v]++;

                    // TODO: LEFTOVER ML 20190702
                    // System.out.println("ins");

                } else if (keyDEL.contains(VCFVarPos)) {
                    tempFinderHM = hmDEL.get(VCFVarPos);
                    indelKS = tempFinderHM.keySet();
                    for (Integer indel : indelKS) {
                        currAltAllele = tempFinderHM.get(indel);
                        if (tempReporterHM.containsKey(currAltAllele)) {
                            currVarObj = tempReporterHM.get(currAltAllele);

                            // TODO: LEFTOVER ML 20190702
                            // VEFFile.print(currVarObj.intM + "=" + currVarObj.varCode + ";");

                            readInfo.append(currVarObj.intM + "=" + currVarObj.varCode + ";");

                            // TODO: LEFTOVER ML 20190702
                            // hmVarCount.put(VCFVarPos, hmVarCount.get(VCFVarPos) + 1);

                            pool_ct[v]++;

                        } else {
                            readInfo.append(VCFVarPos + "=0;");

                            // TODO: LEFTOVER ML 20190702
                            // System.out.println(VCFVarPos + " done");

                        }

                        // TODO: LEFTOVER ML 20190702
                        // System.out.println("del");

                    }
                    total_ct[v]++;
                }
                v++;
            }

            // TODO: LEFTOVER ML 20190702
            // System.out.println(readInfo.toString());

            if (readInfo.toString().isEmpty()) {
                BAMScanner.nextLine(); 	// Skip the QUAL.
                hmMATCH.clear();
                hmINS.clear();
                hmDEL.clear();
                continue;
            }

            // TODO: LEFTOVER ML 20190702
            // System.out.println(readInfo);

            VEFFile.println(QNAME + ":\t" + readInfo + "\t//\t" + startPOS + "\t" + endPOS);
            BAMScanner.nextLine(); 	// Skip the QUAL.
            hmMATCH.clear();
            hmINS.clear();
            hmDEL.clear();
        }

        VEFFile.close();
        BAMScanner.close();
        double[] pool_freq = new double[variantEncoder.size()];
        for (int v = 0; v < numVars; v++) {
            pool_freq[v] = pool_ct[v] / total_ct[v];
        }
        return pool_freq;
    }

}
