import java.util.*;

public class GreedyMotifSearch {

    private static final int NUM_OF_ACIDS = 4;

    public static MotifSearchResponse motifSearch(List<String> dnaList, int numOfSeq, int seqLength, int motifLen) {

        MotifSearchResponse response = new MotifSearchResponse();

        int[] bestMotifSeq = new int[numOfSeq];
        int[] startSeq = new int[numOfSeq];

        /*
        * Initializing sequence arrays
        * */
        for(int i = 0; i < numOfSeq; i++) {
            bestMotifSeq[i] = 0;
            startSeq[i] = 0;
        }

        /*

         // TRADITIONAL WAY
        for(int start1 = 0; start1 < seqLength - motifLen +1; start1++) {
            for(int start2 = 0; start2 < seqLength - motifLen +1; start2++) {
                startSeq[0] = start1;
                startSeq[1] = start2;
                if(MotifSearchUtil.getScore(startSeq, 2, dnaList, motifLen)
                        > MotifSearchUtil.getScore(bestMotifSeq, 2, dnaList, motifLen)) {
                    bestMotifSeq[0] = start1;
                    bestMotifSeq[1] = start2;
                }
            }
        }
        */
        int bestMotifScore = 0;
        int bestMotif1Index = 0;
        int bestMotif2Index = 0;
        int bestMotif1StartVal = 0;
        int bestMotif2StartVal = 0;

        /*
        * Finding first two best motif start index.
        * */

        Random random = new Random();
        random.setSeed(2144567889);
        Set usedSeeds = new HashSet<String>();

        for(int counter = 0; counter < 20; counter++ ){

            for(int i = 0; i < numOfSeq; i++) {
                bestMotifSeq[i] = 0;
                startSeq[i] = 0;
            }

            List seedList = new LinkedList();

            int seedSeq1Index =  Math.abs(random.nextInt(numOfSeq)) % numOfSeq;
            int seedSeq2Index = Math.abs(random.nextInt(numOfSeq)) % numOfSeq;

            while(usedSeeds.contains(seedSeq1Index + "," + seedSeq2Index)
                    || usedSeeds.contains(seedSeq2Index + "," + seedSeq1Index)) {
                seedSeq1Index =  Math.abs(random.nextInt(numOfSeq)) % numOfSeq;
                seedSeq2Index = Math.abs(random.nextInt(numOfSeq)) % numOfSeq;
            }


            usedSeeds.add(seedSeq1Index + "," + seedSeq2Index);
            usedSeeds.add(seedSeq2Index + "," + seedSeq1Index);

            seedList.add(seedSeq1Index);
            seedList.add(seedSeq2Index);


            int bestSeedScore = 0;
            for(int start1 = 0; start1 < seqLength - motifLen +1; start1++) {
                for(int start2 = 0; start2 < seqLength - motifLen +1; start2++) {

                    startSeq[seedSeq1Index] = start1;
                    startSeq[seedSeq2Index] = start2;
                    int seedScore = getScoreForRandomSeeds(startSeq, 2, dnaList, motifLen, seedList);
                    if(seedScore
                            > getScoreForRandomSeeds(bestMotifSeq, 2, dnaList, motifLen, seedList)) {
                        bestMotifSeq[seedSeq1Index] = start1;
                        bestMotifSeq[seedSeq2Index] = start2;
                        bestSeedScore = seedScore;

                    }
                }
            }

            if(bestSeedScore > bestMotifScore) {
                bestMotif1Index = seedSeq1Index;
                bestMotif2Index = seedSeq2Index;
                bestMotif1StartVal = bestMotifSeq[seedSeq1Index];
                bestMotif2StartVal = bestMotifSeq[seedSeq2Index];
                bestMotifScore = bestSeedScore;
            }

        }


        startSeq[bestMotif1Index] = bestMotif1StartVal;
        startSeq[bestMotif2Index] = bestMotif2StartVal;
        bestMotifSeq[bestMotif1Index] = bestMotif1StartVal;
        bestMotifSeq[bestMotif2Index] = bestMotif2StartVal;

        int level = 3;
        List seedSequenceIndex = new LinkedList();
        seedSequenceIndex.add(bestMotif1Index);
        seedSequenceIndex.add(bestMotif2Index);
        for(int index = 0; index < numOfSeq && level <= numOfSeq; index++) {
            if(index == bestMotif1Index || index == bestMotif2Index) {
                continue;
            }
            seedSequenceIndex.add(index);
            for(int sIndex = 0; sIndex < seqLength - motifLen + 1; sIndex++) {
                startSeq[index] = sIndex;
                if(getScoreForRandomSeeds(startSeq, level, dnaList, motifLen, seedSequenceIndex)
                        > getScoreForRandomSeeds(bestMotifSeq, level, dnaList, motifLen, seedSequenceIndex)) {
                    bestMotifSeq[index] = sIndex;
                }
            }
            startSeq[index] = bestMotifSeq[index];
            level++;
        }

        response.setStartSequences(bestMotifSeq);

        /*
        * Calculating the maximum score
        * */
        int maxScore = BBMotifSearch.getScore(bestMotifSeq, numOfSeq, dnaList, motifLen);

        response.setBestMotifScore(maxScore);
        return response;
    }

    private static int getScoreForRandomSeeds(int[] startSeq, int level, List<String> dnaSeq, int motifLen, List<Integer> seedList) {
        int score = 0;
        /*
        * STEP-0: Converting List of String into char[][]
        * */

        char[][] dnaSeqArray = new char[dnaSeq.size()][dnaSeq.get(0).length()];
        int tempIndex = 0;
        for (int seedIndex : seedList) {
            dnaSeqArray[tempIndex] = dnaSeq.get(seedIndex).toCharArray();
            tempIndex++;
        }


        /*
        * STEP1: Building the alignment matrix
        * */

        char[][] alignment = new char[level][motifLen];

        for (int i = 0; i < level; i++) {
            for (int j = 0, dnaCharIndex = startSeq[seedList.get(i)]; j < motifLen
                    && dnaCharIndex < startSeq[seedList.get(i)] + motifLen; j++, dnaCharIndex++) {
                if (dnaCharIndex >= dnaSeq.get(0).length())
                    return 0;
                else
                    alignment[i][j] = dnaSeqArray[i][dnaCharIndex];
            }
        }

        /*
        * STEP2: Building the Profile
        * */
        int[][] profile = new int[NUM_OF_ACIDS][motifLen];

        // Initializing the profile.
        for (int i = 0; i < profile.length; i++)
            for (int j = 0; j < profile[0].length; j++)
                profile[i][j] = 0;


        for (int i = 0; i < alignment.length; i++) {
            for (int j = 0; j < alignment[0].length; j++) {
                int profileIndex = convertDNAAcidCodeToInt(alignment[i][j]) - 1;
                profile[profileIndex][j] = profile[profileIndex][j] + 1;
            }
        }

        /*
        * Calculate Score.
        * */
        int localMax = 0;
        for (int j = 0; j < profile[0].length; j++) {
            localMax = 0;
            for (int i = 0; i < profile.length; i++) {
                localMax = Math.max(localMax, profile[i][j]);
            }
            score += localMax;
        }

        //System.out.println("Score :: " + score);
        return score;

    }
    
    private static int convertDNAAcidCodeToInt(char dnaAcidCode) {
        if (dnaAcidCode == 'A'
                || dnaAcidCode == 'a')
            return 1;
        else if (dnaAcidCode == 'C'
                || dnaAcidCode == 'c')
            return 2;
        else if (dnaAcidCode == 'G'
                || dnaAcidCode == 'g')
            return 3;
        else
            return 4;
    }

}
