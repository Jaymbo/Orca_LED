# Orca_LED

## Installation Instructions

To install the pipeline, you need a computer and a server with "Sbatch" pre-installed and 48 CPUs. The pipeline has been specifically developed for **Justus2**, but it should also work on similar systems.

### Steps

1. **Login to the Server**  
   Use Visual Studio Code (VSC) to log in to the server.  
   For **Justus2**, a port forwarding must be integrated into the login process since the server does not allow external access to webpages. Port forwarding enables local access to the webpage for usage.

2. **Clone the Repository**  
   Navigate to the desired directory and run the following command:  
   ```bash
   git clone https://github.com/Jaymbo/Orca_LED
   ```
3. Set Up the Environment
  Use conda to configure the environment for the pipeline.

4. Start the Web Server
  Log in again, this time directly into the newly created folder from the repository.
  Use Visual Studio Code tasks to automate starting the web server. A task is used to initiate the web server.
