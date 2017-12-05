import org.biojava.nbio.data.sequence.FastaSequence;
import org.biojava.nbio.data.sequence.SequenceUtil;

import java.io.*;
import java.util.List;

public class ProtienReader {


    public static List<FastaSequence> readProtien(String fatsaFileName) throws IOException {

        List<FastaSequence> fatsaSeqList = SequenceUtil.readFasta(new BufferedInputStream(new FileInputStream(fatsaFileName)));

        return fatsaSeqList;
    }
}
