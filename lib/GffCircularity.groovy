import java.nio.file.Files
import java.nio.file.Path

/**
 * Strictly derive mitogenome topology from the Is_circular attribute on GFF
 * region features.
 */
class GffCircularity {

    static Map annotate(Map meta, Collection annotationFiles) {
        def gffFiles = annotationFiles.findAll {
            it.name.toLowerCase().endsWith('.gff')
        }

        if (!gffFiles) {
            throw new IllegalArgumentException(
                "Assembly '${assemblyName(meta)}' has no GFF file; cannot determine Is_circular."
            )
        }
        if (gffFiles.size() > 1) {
            def names = gffFiles.collect { it.name }.sort().join(', ')
            throw new IllegalArgumentException(
                "Assembly '${assemblyName(meta)}' has multiple GFF files (${names}); cannot determine a unique Is_circular value."
            )
        }

        boolean circular = parse(gffFiles.first() as Path, assemblyName(meta))
        return meta + [circular: circular]
    }

    static boolean parse(Path gff, String assembly = null) {
        def source = assembly ?: gff.fileName.toString()
        def verdicts = [] as Set
        int lineNumber = 0

        Files.newBufferedReader(gff).withCloseable { reader ->
            reader.eachLine { line ->
                lineNumber++
                if (!line || line.startsWith('#')) {
                    return
                }

                def fields = line.split('\\t', -1)
                if (fields.size() < 3 || !fields[2].equalsIgnoreCase('region')) {
                    return
                }
                if (fields.size() < 9) {
                    throw invalidValue(source, gff, lineNumber, 'malformed region record')
                }

                def attributes = fields[8].split(';', -1)
                attributes.each { attribute ->
                    def pair = attribute.trim().split('=', 2)
                    if (!pair || !pair[0].equalsIgnoreCase('Is_circular')) {
                        return
                    }
                    if (pair.size() != 2) {
                        throw invalidValue(source, gff, lineNumber, 'missing value')
                    }

                    def value = pair[1].trim()
                    if (value.equalsIgnoreCase('true')) {
                        verdicts << true
                    } else if (value.equalsIgnoreCase('false')) {
                        verdicts << false
                    } else {
                        throw invalidValue(source, gff, lineNumber, "'${value}' is not true or false")
                    }
                }
            }
        }

        if (!verdicts) {
            throw new IllegalArgumentException(
                "Assembly '${source}' GFF '${gff.fileName}' has no valid Is_circular=true|false attribute on a region feature."
            )
        }
        if (verdicts.size() > 1) {
            throw new IllegalArgumentException(
                "Assembly '${source}' GFF '${gff.fileName}' has conflicting Is_circular values on region features."
            )
        }
        return verdicts.first() as boolean
    }

    private static IllegalArgumentException invalidValue(String assembly, Path gff, int line, String reason) {
        return new IllegalArgumentException(
            "Assembly '${assembly}' GFF '${gff.fileName}' has invalid Is_circular at line ${line}: ${reason}."
        )
    }

    private static String assemblyName(Map meta) {
        return (meta.mt_assembly_prefix ?: meta.id ?: 'unknown').toString()
    }
}
