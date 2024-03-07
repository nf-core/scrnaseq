import groovy.json.JsonSlurper


class WorkflowScrnaseq {
    // Retrieve the aligner-specific protocol based on the specified protocol.
    // Returns a map ["protocol": protocol, "extra_args": <extra args>, "whitelist": <path to whitelist>]
    // extra_args and whitelist are optional.
    public static Map getProtocol(workflow, log, aligner, protocol) {
        def jsonSlurper = new JsonSlurper()
        def json = new File("${workflow.projectDir}/assets/protocols.json").text
        def protocols = jsonSlurper.parseText(json)
        def aligner_map = protocols[aligner]
        if(aligner_map.containsKey(protocol)) {
            return aligner_map[protocol]
        } else {
            log.warn("Protocol '${protocol}' not recognized by the pipeline. Passing on the protocol to the aligner unmodified.")
            return ["protocol": protocol]
        }
    }

}
