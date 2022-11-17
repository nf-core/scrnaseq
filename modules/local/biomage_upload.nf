process BIOMAGE_UPLOAD {
    label 'process_low'
    container 'biomage/programmatic-interface:0.0.7'

    input:
    val email
    val password
    val instance_url
    path samples

    output:
    stdout

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python
    import biomage_programmatic_interface as bpi

    connection = bpi.Connection('$email', '$password', '$instance_url', verbose=False)
    experiment_id = connection.create_experiment()
    for sample in '$samples'.split():
        connection.upload_samples(experiment_id, sample)

    print(f"Project {experiment_id} successfuly created!")
    print(f'You can view it at https://$instance_url/')
    """
}
