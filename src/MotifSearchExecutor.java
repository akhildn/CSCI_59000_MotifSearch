import org.biojava.nbio.data.sequence.FastaSequence;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class MotifSearchExecutor {



    public static void main(String args[]) {

        System.out.println("Hello World");

        try {

            /*
            * Reading DNA sequences as List of Strings
            * */
            String fastaFileName = args[0];
            int motifLength = Integer.parseInt(args[1]);
            String algoType = args[2];


            List<String> dnaList = getDNASeqStrings(fastaFileName);

            if(algoType.equalsIgnoreCase("bb")) {
                System.out.println("Branch Bound ::");
                invokeAndMeasureBBMotifSearch(dnaList, motifLength);
            }
            else if(algoType.equalsIgnoreCase("greedy")) {
                System.out.println("Modified Greedy ::");
                invokeAndMeasureGreedyMotifSearch(dnaList, motifLength);
            } else {
                System.out.println("Not a Valid Algorithm name passed.");
            }



        } catch (IOException e) {
            System.out.println("Error reading the Fasta File");
            e.printStackTrace();
        } catch (NumberFormatException ex) {
            System.out.println("Motif Length should be integer");
            ex.printStackTrace();
        } catch(Exception e) {
            e.printStackTrace();
        }

    }


    private static void invokeAndMeasureGreedyMotifSearch(List<String> dnaList, int motifLength) {
        MotifSearchResponse response;
        double timeTaken;
        long startGreedyTime = System.currentTimeMillis();
        response = GreedyMotifSearch.motifSearch(dnaList, dnaList.size(), dnaList.get(0).length(), motifLength);
        System.out.print("Start Sequences :: [");
        for (int startIndex : response.getStartSequences()) {
            System.out.print(startIndex + ",");
        }
        System.out.println("] \t Max Score = " + response.getBestMotifScore());
        long endGreedyTime = System.currentTimeMillis();
        timeTaken = (endGreedyTime - startGreedyTime);

        System.out.println("\n Time Taken for Greedy: " + timeTaken);
    }

    private static void invokeAndMeasureBBMotifSearch(List<String> dnaList, int motifLength) {
        long startBBTime = System.currentTimeMillis();
        MotifSearchResponse response = BBMotifSearch.branchBoundMotifSearch(dnaList, dnaList.size()
                                                                        , dnaList.get(0).length(), motifLength);
        System.out.print("Start Sequences :: [");
        for (int startIndex : response.getStartSequences()) {
            System.out.print(startIndex + ",");
        }
        System.out.println("] \t Max Score = " + response.getBestMotifScore());
        long endBBTime = System.currentTimeMillis();
        double timeTaken = (endBBTime - startBBTime);

        System.out.println("\n Time Taken for BB: " + timeTaken);
    }

    private static List<String> getDNASeqStrings(String fastaFileName) throws IOException {
        List<FastaSequence> fastaSeq = ProtienReader.readProtien(fastaFileName);
        List<String> dnaList = new LinkedList<>();
        for (FastaSequence sequence : fastaSeq) {
            dnaList.add(sequence.getSequence());
        }
        return dnaList;
    }

}
