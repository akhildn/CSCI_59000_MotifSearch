import org.biojava.nbio.data.sequence.FastaSequence;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class MotifSearch {

    private static final int NUM_OF_ACIDS = 4;
    public static void main(String args[]) {

        System.out.println("Hello World");

        try {
            List<FastaSequence> fastaSeq = ProtienReader.readProtien("test0.fasta");
            List<String> dnaList = new ArrayList<>();
            for(FastaSequence sequence : fastaSeq) {
               dnaList.add(sequence.getSequence());
            }


            System.out.println("Finished Reading");

          /*  dnaList.add("GGCGTTCAGGCA");
            dnaList.add("AAGAATCAGTCA");
            dnaList.add("CAAGGAGTTCGC");
            dnaList.add("CACGTCAATCAC");
            dnaList.add("CAATAATATTCG");*/

            MotifSearchResponse response = branchBoundMotifSearch(dnaList, dnaList.size(), dnaList.get(0).length(), 8);

            System.out.println("Start Sequences ::" + response.getLevel());
            for(int startIndex : response.getStartSequences()) {
                System.out.print(startIndex + "\t");
            }

        } catch (IOException e) {
            System.out.println("Error reading the Fasta File");
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

    }

    public static MotifSearchResponse branchBoundMotifSearch(List<String> dnaSeq, int numOfSeq, int lengthOfSeq, int motifLen) throws InterruptedException {

        int[] startSequences = new int[numOfSeq];
        int[] bestSequences = new int[numOfSeq];
        for(int i = 0; i < startSequences.length; i++) {
            startSequences[i] = 0;
        }

        MotifSearchResponse response;
        MotifSearchResponse bestMotifResponse = new MotifSearchResponse();
        int bestScore = 0;
        int level = 1;

        while(level > 0) {
            if(level < numOfSeq) {
                //System.out.println("Level :: " + level);
                int optimisticScore = getScore(startSequences, level, dnaSeq, motifLen) + (numOfSeq - level) * motifLen;
                if(optimisticScore < bestScore) {
                    response = byPass(startSequences, level, numOfSeq, lengthOfSeq - motifLen + 1);
                    startSequences = response.getStartSequences();
                    level = response.getLevel();
                }
                else {
                    response = nextVertex(startSequences, level, numOfSeq, lengthOfSeq - motifLen + 1);
                    startSequences = response.getStartSequences();
                    level = response.getLevel();
                }
            } else {
                //System.out.println("Inside else :: " + level);
                int score = getScore(startSequences, numOfSeq, dnaSeq, motifLen);
                if(score > bestScore) {
                    bestScore = score;
                    for(int i = 0; i < bestSequences.length; i++) {
                        bestSequences[i] = startSequences[i];
                    }
                    bestMotifResponse.setStartSequences(bestSequences);
                    bestMotifResponse.setBestMotifScore(bestScore);
                    bestMotifResponse.setLevel(level);
                }

                response = nextVertex(startSequences, level, numOfSeq, lengthOfSeq - motifLen + 1);
                //System.out.println("After response :: " + response.getLevel());
                startSequences = response.getStartSequences();
                level = response.getLevel();
            }
        }

        return bestMotifResponse;
    }


    private static MotifSearchResponse nextVertex(int[] startSequences, int level, int numOfSeq, int i) throws InterruptedException {
        //Thread.sleep(10);
        MotifSearchResponse response = new MotifSearchResponse();
        if (level < numOfSeq) {
            //System.out.println("IF :: " + level);
            startSequences[level + 1 -1] = 0;
            response.setStartSequences(startSequences);
            response.setLevel(level + 1);
            return response;
        } else {
           // System.out.println(startSequences + "::" + level + "::" + numOfSeq + "::" + i);
            for (int j = numOfSeq; j > 0; j--) {
                if (startSequences[j-1] < i-1) {
                   // System.out.println(startSequences[j-1] + "::" + level + "::" + numOfSeq + "::" + (i-1));
                    startSequences[j-1] += 1;
                    response.setLevel(j);
                    response.setStartSequences(startSequences);

                    return response;
                }
            }
        }
        response.setLevel(0);
        response.setStartSequences(startSequences);
        return response;

    }

    private static MotifSearchResponse byPass(int[] startSequences, int level, int numOfSeq, int k) {
        //System.out.println("ByPass");
        MotifSearchResponse response = new MotifSearchResponse();
        for (int i = level; i > 0; i--) {
            if (startSequences[i-1] < k-1) {
                startSequences[i-1] += 1;
                response.setStartSequences(startSequences);
                response.setLevel(i);
                return response;
            }
        }
        response.setStartSequences(startSequences);
        response.setLevel(0);
        return response;
    }

    private static int getScore(int[] startSeq, int level, List<String> dnaSeq, int motifLen) {

        int score = 0;
        /*
        * STEP-0: Converting List of String into char[][]
        * */

        char[][] dnaSeqArray = new char[dnaSeq.size()][dnaSeq.get(0).length()];
        int tempIndex = 0;
        for(String seq : dnaSeq) {
            dnaSeqArray[tempIndex] = dnaSeq.get(tempIndex).toCharArray();
            tempIndex++;
        }


        /*
        * STEP1: Building the alignment matrix
        * */

        char[][] alignment = new char[level][motifLen];

        for(int i = 0; i < level; i++) {
            for(int j = 0, dnaCharIndex = startSeq[i]; j < motifLen
                    && dnaCharIndex < startSeq[i] + motifLen; j++, dnaCharIndex++) {
                if(dnaCharIndex >= dnaSeq.get(0).length())
                    return 0;
                else
                    alignment[i][j] = dnaSeqArray[i][dnaCharIndex];
            }
        }

/*        System.out.println("DEBUG");

        for(int i = 0; i < alignment.length; i++) {
            for(int j = 0; j < alignment[i].length; j++) {
                System.out.print(alignment[i][j] + "\t");
            }
            System.out.println();
        }*/


        /*
        * STEP2: Building the Profile
        * */
        int[][] profile = new int[NUM_OF_ACIDS][motifLen];

        // Initializing the profile.
        for(int i = 0; i < profile.length; i++)
            for(int j = 0; j < profile[0].length; j++)
                profile[i][j] = 0;


        for(int i = 0; i < alignment.length; i++) {
            for(int j = 0; j < alignment[0].length; j++) {
                int profileIndex = convertDNAAcidCodeToInt(alignment[i][j]) - 1;
                profile[profileIndex][j] = profile[profileIndex][j] + 1;
            }
        }

/*        System.out.println("DEBUG");

        for(int i = 0; i < profile.length; i++) {
            for(int j = 0; j < profile[i].length; j++) {
                System.out.print(profile[i][j] + "\t");
            }
            System.out.println();
        }*/

        /*
        * Calculate Score.
        * */
        int localMax = 0;
        for(int j = 0; j < profile[0].length; j++) {
            localMax = 0;
            for(int i = 0; i < profile.length; i++) {
                localMax = Math.max(localMax, profile[i][j]);
            }
            score += localMax;
        }

        //System.out.println("Score :: " + score);
        return score;
    }

    private static int convertDNAAcidCodeToInt(char dnaAcidCode) {
        if(dnaAcidCode == 'A'
                || dnaAcidCode == 'a')
            return 1;
        else if(dnaAcidCode == 'C'
                || dnaAcidCode == 'c')
            return 2;
        else if(dnaAcidCode == 'G'
                || dnaAcidCode == 'g')
            return 3;
        else
            return 4;
    }
}
