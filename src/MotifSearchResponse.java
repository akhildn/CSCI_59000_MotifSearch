/**
 * Class to represent the best motif and its score.
 */
public class MotifSearchResponse {
    private int[] startSequences;
    private int bestMotifScore;
    private int level;

    /**
     * @return the best motif score
     */
    public int getBestMotifScore() {
        return bestMotifScore;
    }

    /**
     * @param bestMotifScore
     */
    public void setBestMotifScore(int bestMotifScore) {
        this.bestMotifScore = bestMotifScore;
    }

    /**
     * @return the best motif sequence
     */
    public int[] getStartSequences() {
        return startSequences;
    }

    /**
     * @param startSequences
     */
    public void setStartSequences(int[] startSequences) {
        this.startSequences = startSequences;
    }

    public int getLevel() {
        return level;
    }

    public void setLevel(int level) {
        this.level = level;
    }
}
