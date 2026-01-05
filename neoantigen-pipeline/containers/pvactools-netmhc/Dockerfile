# Dockerfile for pVACtools with NetMHCpan/NetMHCIIpan
# Base: griffithlab/pvactools with custom NetMHC installations

FROM griffithlab/pvactools:4.1.1

LABEL maintainer="Bioinformatics Core Facility"
LABEL description="pVACtools with NetMHCpan and NetMHCIIpan"

# Install additional dependencies
RUN apt-get update && apt-get install -y \
    tcsh \
    gawk \
    && rm -rf /var/lib/apt/lists/*

# Create directory for NetMHC tools
RUN mkdir -p /opt/netmhc

# NetMHCpan and NetMHCIIpan must be downloaded separately with academic license
# The archives should be placed in the build context directory

# Copy and extract NetMHCpan (you need to download this separately)
# COPY netMHCpan-4.1.Linux.tar.gz /opt/netmhc/
# RUN cd /opt/netmhc && tar -xzf netMHCpan-4.1.Linux.tar.gz && rm netMHCpan-4.1.Linux.tar.gz

# Copy and extract NetMHCIIpan (you need to download this separately)  
# COPY netMHCIIpan-4.0.Linux.tar.gz /opt/netmhc/
# RUN cd /opt/netmhc && tar -xzf netMHCIIpan-4.0.Linux.tar.gz && rm netMHCIIpan-4.0.Linux.tar.gz

# Note: After uncommenting the above and building, you'll need to:
# 1. Edit /opt/netmhc/netMHCpan-4.1/netMHCpan to set NMHOME
# 2. Edit /opt/netmhc/netMHCIIpan-4.0/netMHCIIpan to set NMHOME
# 3. Download and extract the data files

# Placeholder script that shows instructions
RUN echo '#!/bin/bash\n\
echo "NetMHCpan/NetMHCIIpan not installed"\n\
echo "Please build this container with your licensed NetMHC tools"\n\
echo "See Dockerfile for instructions"\n\
exit 1' > /opt/netmhc/placeholder.sh && chmod +x /opt/netmhc/placeholder.sh

# Environment variables (update paths after installation)
ENV PATH="/opt/netmhc/netMHCpan-4.1:/opt/netmhc/netMHCIIpan-4.0:${PATH}"
ENV TMPDIR="/tmp"

# IEDB directory for pVACtools
RUN mkdir -p /opt/iedb

# Set working directory
WORKDIR /data

# Default command
CMD ["pvacseq", "--help"]
