import org.biojava.nbio.data.sequence.FastaSequence;

import java.io.IOException;
import java.util.List;

public class MotifSearch {


    public static void main(String args[]) {

        System.out.println("Hello World");

        try {
            List<FastaSequence> fastaSeq = ProtienReader.readProtien("Sample.fasta");

            for(FastaSequence sequence : fastaSeq) {
                System.out.println("New Sequence \n" + sequence.getId() + "\n" + sequence.getSequence());
            }


            System.out.println("Finished Reading");
        } catch (IOException e) {
            System.out.println("Error reading the Fasta File");
            e.printStackTrace();
        }

    }
}
