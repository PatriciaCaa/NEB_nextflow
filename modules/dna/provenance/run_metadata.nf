process CAPTURE_RUN_METADATA {

    tag "DNA reproducibility metadata"

    label 'python'

    publishDir "${params.results_dir}/reproducibility", mode: 'copy'

    input:
    val metadata_json

    output:
    path "run_metadata.json", emit: metadata

    script:
    // Encode the metadata to preserve JSON content when passing it to the shell.

    def encoded = metadata_json.bytes.encodeBase64().toString()
    """
    echo '${encoded}' | base64 -d > run_metadata.json

    python3 -m json.tool run_metadata.json > /dev/null
    """
}
