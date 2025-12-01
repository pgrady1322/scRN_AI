# Phase 3 Implementation Plan: Configuration & Orchestration

**Status**: 🚧 IN PROGRESS  
**Start Date**: November 25, 2025  
**Estimated Duration**: 1-2 weeks  
**Objective**: Transform individual CLI commands into a fully automated, config-driven pipeline system

---

## 📋 Overview

Phase 3 adds automation and orchestration capabilities to scrn_ai, enabling:
- YAML-driven configuration instead of manual command parameters
- Automatic workflow chaining (preprocess → normalize → umap → pseudotime)
- Docker Compose multi-service orchestration
- Pipeline checkpointing and resume capabilities
- Intermediate file management

---

## 🎯 Implementation Roadmap

### Milestone 1: YAML Configuration Parser (Days 1-3)
**Deliverable**: `scrn_ai/config/parser.py`

#### Tasks:
- [x] Create `scrn_ai/config/` directory
- [ ] Implement `ConfigParser` class
- [ ] Add YAML validation (schema checking)
- [ ] Support command-line overrides
- [ ] Add environment variable substitution
- [ ] Create default configuration template
- [ ] Add error handling and informative messages
- [ ] Write unit tests for parser

#### Components to Build:

**File**: `scrn_ai/config/parser.py`
```python
class ConfigParser:
    """Parse and validate config.yaml files"""
    
    def __init__(self, config_path: str):
        """Load configuration from YAML file"""
        
    def validate(self) -> bool:
        """Validate configuration against schema"""
        
    def get_input_config(self) -> dict:
        """Extract input configuration"""
        
    def get_preprocessing_config(self) -> dict:
        """Extract preprocessing parameters"""
        
    def get_normalization_config(self) -> dict:
        """Extract normalization parameters"""
        
    def get_analysis_config(self) -> dict:
        """Extract analysis configuration"""
        
    def merge_cli_overrides(self, cli_args: dict):
        """Override config values with CLI arguments"""
        
    def substitute_env_vars(self):
        """Replace ${VAR} with environment variables"""
```

**File**: `scrn_ai/config/schema.yaml`
```yaml
# Configuration schema for validation
required_fields:
  - input.matrix_path
  - output.results_dir

optional_fields:
  input:
    - metadata_path
    - input_format
  preprocessing:
    - min_genes_per_cell
    - min_cells_per_gene
    # ... etc
```

**File**: `scrn_ai/config/defaults.yaml`
```yaml
# Default configuration values
preprocessing:
  min_genes_per_cell: 200
  min_cells_per_gene: 3
  max_genes_per_cell: null
  max_mito_pct: null

normalization:
  method: "seurat"
  algorithm: "LogNormalize"
  scale_factor: 10000

analysis:
  run_umap: true
  umap_n_neighbors: 15
  umap_min_dist: 0.1
  run_pseudotime: false
  # ... etc
```

#### Testing:
```bash
# Test config parser
python -m pytest tests/test_config_parser.py

# Test with sample config
scrn_ai --config examples/sample_config.yaml --validate-only
```

---

### Milestone 2: Workflow Orchestrator (Days 4-6)
**Deliverable**: Enhanced `scrn_ai/main.py` with pipeline execution

#### Tasks:
- [ ] Create `WorkflowOrchestrator` class
- [ ] Implement workflow DAG (Directed Acyclic Graph)
- [ ] Add automatic step dependencies
- [ ] Implement intermediate file management
- [ ] Add progress tracking and logging
- [ ] Error handling with rollback support
- [ ] Add dry-run mode (show what would execute)
- [ ] Write integration tests

#### Components to Build:

**File**: `scrn_ai/orchestrator.py`
```python
class WorkflowOrchestrator:
    """Orchestrate multi-step analysis pipeline"""
    
    def __init__(self, config: ConfigParser):
        """Initialize with parsed configuration"""
        
    def build_dag(self) -> WorkflowDAG:
        """Build workflow dependency graph"""
        
    def execute(self, dry_run: bool = False):
        """Execute complete pipeline"""
        
    def execute_step(self, step: str):
        """Execute single pipeline step"""
        
    def check_dependencies(self, step: str) -> bool:
        """Verify all dependencies are satisfied"""
        
    def get_execution_plan(self) -> List[str]:
        """Return ordered list of steps to execute"""
```

**File**: `scrn_ai/workflows/dag.py`
```python
class WorkflowDAG:
    """Workflow dependency graph"""
    
    DEPENDENCIES = {
        'preprocess': [],
        'normalize': ['preprocess'],
        'aitype_pre': ['normalize'],  # Optional
        'umap': ['normalize'],
        'pseudotime': ['normalize'],
        'aitype_post': ['pseudotime'],  # Optional
    }
    
    def topological_sort(self) -> List[str]:
        """Return execution order"""
```

**Enhanced**: `scrn_ai/main.py`
```python
@click.command()
@click.option('--config', '-c', type=click.Path(exists=True),
              help='Path to config.yaml file')
@click.option('--dry-run', is_flag=True,
              help='Show execution plan without running')
@click.option('--resume', is_flag=True,
              help='Resume from last checkpoint')
def run_pipeline(config, dry_run, resume):
    """Execute complete analysis pipeline from config file"""
    
    # Parse configuration
    parser = ConfigParser(config)
    
    # Create orchestrator
    orchestrator = WorkflowOrchestrator(parser)
    
    # Execute pipeline
    if dry_run:
        plan = orchestrator.get_execution_plan()
        click.echo("Execution plan:")
        for step in plan:
            click.echo(f"  {step}")
    else:
        orchestrator.execute(resume=resume)
```

#### Usage:
```bash
# Run full pipeline from config
scrn_ai --config config.yaml

# Dry run (show plan)
scrn_ai --config config.yaml --dry-run

# Resume from checkpoint
scrn_ai --config config.yaml --resume
```

---

### Milestone 3: Docker Compose Integration (Days 7-8)
**Deliverable**: Working `docker-compose.yml` with multi-service setup

#### Tasks:
- [ ] Create `docker-compose.yml` template
- [ ] Define service dependencies
- [ ] Set up volume mappings
- [ ] Configure environment variables
- [ ] Add health checks
- [ ] Create example configurations
- [ ] Document docker-compose usage
- [ ] Test multi-stage pipeline execution

#### Files to Create:

**File**: `docker-compose.yml`
```yaml
version: "3.8"

services:
  # Full pipeline with config file
  pipeline:
    build: .
    image: scrn_ai:0.1
    container_name: scrn_ai_pipeline
    volumes:
      - ./data:/data
      - ./config:/config
      - ./output:/output
    environment:
      - OPENAI_API_KEY=${OPENAI_API_KEY}
    command: ["--config", "/config/config.yaml"]
    
  # Individual step: Preprocessing
  preprocess:
    image: scrn_ai:0.1
    container_name: scrn_ai_preprocess
    volumes:
      - ./data:/data
      - ./output:/output
    command:
      - "preprocess"
      - "--input"
      - "/data/input/dataset.h5ad"
      - "--output"
      - "/output/processed.h5ad"
      - "--min-genes"
      - "200"
      - "--min-cells"
      - "3"
    
  # Individual step: Normalization
  normalize:
    image: scrn_ai:0.1
    container_name: scrn_ai_normalize
    depends_on:
      - preprocess
    volumes:
      - ./output:/output
    command:
      - "normalize"
      - "--input"
      - "/output/processed.h5ad"
      - "--output"
      - "/output/normalized.h5ad"
      - "--method"
      - "seurat"
      - "--algorithm"
      - "LogNormalize"
  
  # Individual step: AI Cell Typing (optional)
  aitype:
    image: scrn_ai:0.1
    container_name: scrn_ai_aitype
    depends_on:
      - normalize
    volumes:
      - ./output:/output
    environment:
      - OPENAI_API_KEY=${OPENAI_API_KEY}
    command:
      - "aitype"
      - "--input"
      - "/output/normalized.h5ad"
      - "--output"
      - "/output/cell_types/"
      - "--timing"
      - "pre_analysis"
      - "--model"
      - "gpt-4"
    profiles:
      - ai  # Only run when 'ai' profile is active
  
  # Individual step: UMAP
  umap:
    image: scrn_ai:0.1
    container_name: scrn_ai_umap
    depends_on:
      - normalize
    volumes:
      - ./output:/output
    command:
      - "umap"
      - "--input"
      - "/output/normalized.h5ad"
      - "--output"
      - "/output/umap.png"
      - "--color-by"
      - "leiden"
  
  # Individual step: Pseudotime
  pseudotime:
    image: scrn_ai:0.1
    container_name: scrn_ai_pseudotime
    depends_on:
      - normalize
    volumes:
      - ./output:/output
    command:
      - "pseudotime"
      - "--input"
      - "/output/normalized.h5ad"
      - "--output"
      - "/output/pseudotime/"
      - "--method"
      - "dpt"
      - "--scale"
      - "small"

volumes:
  data:
  output:
  config:
```

**File**: `docker-compose.override.yml` (for local development)
```yaml
version: "3.8"

services:
  pipeline:
    build:
      context: .
      dockerfile: Dockerfile
    volumes:
      - ./scrn_ai:/opt/scrn_ai/scrn_ai:ro  # Mount source for development
```

**File**: `.env.example`
```bash
# OpenAI API Configuration
OPENAI_API_KEY=your-api-key-here

# Data paths
DATA_DIR=./data
OUTPUT_DIR=./output
CONFIG_DIR=./config

# Pipeline settings
SCRN_AI_LOG_LEVEL=INFO
SCRN_AI_CHECKPOINT_DIR=./checkpoints
```

#### Usage:
```bash
# Run full pipeline
docker-compose up pipeline

# Run individual steps in sequence
docker-compose up preprocess normalize umap pseudotime

# Run with AI typing
docker-compose --profile ai up

# Run in detached mode
docker-compose up -d

# View logs
docker-compose logs -f pipeline

# Clean up
docker-compose down -v
```

---

### Milestone 4: Pipeline State Management (Days 9-10)
**Deliverable**: Checkpointing and resume capabilities

#### Tasks:
- [ ] Create `StateManager` class
- [ ] Implement checkpoint creation
- [ ] Add resume from checkpoint logic
- [ ] Track completed steps
- [ ] Store intermediate file paths
- [ ] Add checkpoint cleanup
- [ ] Implement step skipping for completed work
- [ ] Write tests for state management

#### Components to Build:

**File**: `scrn_ai/state.py`
```python
class StateManager:
    """Manage pipeline execution state and checkpoints"""
    
    def __init__(self, checkpoint_dir: str):
        """Initialize state manager"""
        
    def save_checkpoint(self, step: str, outputs: dict):
        """Save checkpoint after successful step completion"""
        
    def load_checkpoint(self) -> dict:
        """Load most recent checkpoint"""
        
    def get_completed_steps(self) -> List[str]:
        """Return list of completed steps"""
        
    def should_skip_step(self, step: str) -> bool:
        """Check if step can be skipped (already completed)"""
        
    def get_intermediate_file(self, step: str, key: str) -> str:
        """Retrieve path to intermediate file from step"""
        
    def clear_checkpoints(self):
        """Remove all checkpoint files"""
```

**File Format**: `.scrn_ai_checkpoint.json`
```json
{
  "version": "0.1.0",
  "timestamp": "2025-11-25T10:30:00Z",
  "config_hash": "abc123...",
  "completed_steps": [
    {
      "step": "preprocess",
      "timestamp": "2025-11-25T10:25:00Z",
      "outputs": {
        "processed_file": "/output/processed.h5ad"
      },
      "duration_seconds": 45.2
    },
    {
      "step": "normalize",
      "timestamp": "2025-11-25T10:28:00Z",
      "outputs": {
        "normalized_file": "/output/normalized.h5ad"
      },
      "duration_seconds": 120.5
    }
  ],
  "current_step": "umap",
  "total_duration_seconds": 165.7
}
```

#### Integration with Orchestrator:
```python
class WorkflowOrchestrator:
    def __init__(self, config: ConfigParser, resume: bool = False):
        self.state = StateManager(config.get('checkpoint_dir', './checkpoints'))
        
        if resume:
            checkpoint = self.state.load_checkpoint()
            self.completed_steps = checkpoint['completed_steps']
    
    def execute_step(self, step: str):
        # Skip if already completed
        if self.state.should_skip_step(step):
            click.echo(f"✓ Skipping {step} (already completed)")
            return
        
        # Execute step
        outputs = self._run_step(step)
        
        # Save checkpoint
        self.state.save_checkpoint(step, outputs)
```

---

## 📁 Directory Structure After Phase 3

```
scRN_AI/
├── scrn_ai/
│   ├── config/                    # NEW
│   │   ├── __init__.py
│   │   ├── parser.py              # Configuration parser
│   │   ├── schema.yaml            # Validation schema
│   │   └── defaults.yaml          # Default values
│   ├── orchestrator.py            # NEW - Workflow orchestration
│   ├── state.py                   # NEW - State management
│   ├── workflows/
│   │   ├── dag.py                 # NEW - Dependency graph
│   │   └── ... (existing workflows)
│   └── ... (existing files)
├── docker-compose.yml             # NEW
├── docker-compose.override.yml    # NEW
├── .env.example                   # NEW
├── examples/                      # NEW
│   ├── sample_config.yaml
│   └── advanced_config.yaml
├── tests/                         # NEW
│   ├── test_config_parser.py
│   ├── test_orchestrator.py
│   └── test_state_manager.py
└── ... (existing files)
```

---

## 🧪 Testing Strategy

### Unit Tests
```bash
# Test configuration parser
pytest tests/test_config_parser.py

# Test orchestrator
pytest tests/test_orchestrator.py

# Test state manager
pytest tests/test_state_manager.py
```

### Integration Tests
```bash
# Test full pipeline with sample data
scrn_ai --config examples/sample_config.yaml --dry-run

# Test resume capability
scrn_ai --config examples/sample_config.yaml
# Kill mid-execution, then:
scrn_ai --config examples/sample_config.yaml --resume
```

### Docker Compose Tests
```bash
# Build images
docker-compose build

# Test individual services
docker-compose up preprocess
docker-compose up normalize

# Test full pipeline
docker-compose up pipeline
```

---

## 📝 Documentation Updates Required

### Files to Update:
1. **README.md**:
   - Add "Pipeline Mode" section
   - Document `--config` flag usage
   - Add docker-compose examples
   - Document checkpoint/resume functionality

2. **PHASE3_COMPLETE.md**:
   - Create completion summary
   - Document all new features
   - Add usage examples
   - Include troubleshooting guide

3. **DEVELOPMENT_ROADMAP.md**:
   - Mark Phase 3 as complete
   - Update Phase 4 planning

---

## 🎯 Success Criteria

Phase 3 is complete when:

- [ ] ✅ Config parser loads and validates YAML files
- [ ] ✅ Command-line overrides work correctly
- [ ] ✅ Environment variable substitution functions
- [ ] ✅ Orchestrator executes full pipeline from config
- [ ] ✅ Workflow DAG correctly orders steps
- [ ] ✅ Docker Compose runs multi-service pipeline
- [ ] ✅ Checkpointing saves intermediate state
- [ ] ✅ Resume successfully continues from checkpoint
- [ ] ✅ All tests pass (unit + integration)
- [ ] ✅ Documentation is complete and accurate
- [ ] ✅ Example configurations work end-to-end

---

## 📊 Timeline

| Days | Milestone | Deliverables |
|------|-----------|--------------|
| 1-3 | Config Parser | `parser.py`, `schema.yaml`, `defaults.yaml`, tests |
| 4-6 | Orchestrator | `orchestrator.py`, `dag.py`, enhanced `main.py` |
| 7-8 | Docker Compose | `docker-compose.yml`, `.env.example`, docs |
| 9-10 | State Management | `state.py`, checkpoint system, resume logic |
| 11-12 | Testing & Docs | Integration tests, documentation, examples |
| 13-14 | Buffer/Polish | Bug fixes, optimization, final testing |

**Total**: ~2 weeks (14 days)

---

## 🚀 Getting Started

### Step 1: Set Up Directory Structure
```bash
mkdir -p scrn_ai/config
mkdir -p examples
mkdir -p tests
touch scrn_ai/config/__init__.py
```

### Step 2: Create Configuration Schema
Start with `scrn_ai/config/schema.yaml` based on README examples.

### Step 3: Implement Config Parser
Build `scrn_ai/config/parser.py` with validation logic.

### Step 4: Test Config Parser
Create sample configs and test parsing.

### Continue through milestones...

---

## 📌 Notes

- Maintain backward compatibility with existing CLI usage
- Config file should be optional (CLI-only mode still works)
- Consider performance for large datasets (lazy loading)
- Add comprehensive logging for debugging
- Document all new environment variables
- Provide migration guide for existing users

---

**Ready to start? Let's begin with Milestone 1: Configuration Parser! 🎉**
