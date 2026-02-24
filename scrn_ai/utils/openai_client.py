#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRN_AI v0.1.0

OpenAI API wrapper for cell type prediction with rate limiting.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: GNU General Public License v3.0 - See LICENSE
"""

import os
import time
import json
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass

try:
    from openai import OpenAI
    OPENAI_AVAILABLE = True
except ImportError:
    OPENAI_AVAILABLE = False
    print("Warning: openai package not installed. Install with: pip install openai")


@dataclass
class CellTypePrediction:
    """Container for cell type prediction results."""
    cluster_id: str
    predicted_type: str
    confidence: float
    reasoning: str
    marker_genes_used: List[str]
    alternative_types: List[Tuple[str, float]] = None  # [(type, confidence), ...]


class OpenAIClient:
    """
    Client for OpenAI API interactions focused on cell type identification.
    
    Handles:
    - API authentication
    - Rate limiting
    - Prompt construction
    - Response parsing
    - Error handling
    """
    
    def __init__(
        self,
        api_key: Optional[str] = None,
        model: str = "gpt-4",
        temperature: float = 0.3,
        max_retries: int = 3,
        retry_delay: int = 2
    ):
        """
        Initialize OpenAI client.
        
        Parameters
        ----------
        api_key : str, optional
            OpenAI API key. If None, reads from OPENAI_API_KEY env variable.
        model : str
            Model to use (gpt-4, gpt-4-turbo, gpt-3.5-turbo).
        temperature : float
            Sampling temperature (0.0-1.0). Lower = more deterministic.
        max_retries : int
            Maximum number of API call retries on failure.
        retry_delay : int
            Seconds to wait between retries.
        """
        if not OPENAI_AVAILABLE:
            raise ImportError("openai package not installed. Install with: pip install openai")
        
        # Get API key from parameter or environment
        self.api_key = api_key or os.getenv("OPENAI_API_KEY")
        if not self.api_key:
            raise ValueError(
                "OpenAI API key not found. Set OPENAI_API_KEY environment variable "
                "or pass api_key parameter."
            )
        
        self.client = OpenAI(api_key=self.api_key)
        self.model = model
        self.temperature = temperature
        self.max_retries = max_retries
        self.retry_delay = retry_delay
        
        # Rate limiting (adjust based on your API tier)
        self.last_call_time = 0
        self.min_call_interval = 1.0  # seconds between calls
    
    def _rate_limit(self):
        """Enforce rate limiting between API calls."""
        current_time = time.time()
        time_since_last_call = current_time - self.last_call_time
        
        if time_since_last_call < self.min_call_interval:
            sleep_time = self.min_call_interval - time_since_last_call
            time.sleep(sleep_time)
        
        self.last_call_time = time.time()
    
    def _build_cell_type_prompt(
        self,
        marker_genes: List[str],
        cluster_id: str,
        species: str = "human",
        tissue: Optional[str] = None,
        context: Optional[str] = None
    ) -> str:
        """
        Build a prompt for cell type identification.
        
        Parameters
        ----------
        marker_genes : list of str
            Top marker genes for the cluster.
        cluster_id : str
            Cluster identifier.
        species : str
            Species (human, mouse, etc.).
        tissue : str, optional
            Tissue type if known.
        context : str, optional
            Additional context (e.g., developmental stage, disease state).
        
        Returns
        -------
        str
            Formatted prompt for the LLM.
        """
        tissue_info = f" from {tissue} tissue" if tissue else ""
        context_info = f"\n\nAdditional context: {context}" if context else ""
        
        prompt = f"""You are an expert in single-cell RNA sequencing data analysis and cell type identification.

I have a cluster (Cluster {cluster_id}) from a {species} single-cell RNA-seq dataset{tissue_info} with the following top marker genes:

{', '.join(marker_genes)}

Based on these marker genes, please:

1. Identify the most likely cell type for this cluster
2. Provide your confidence level (0.0-1.0)
3. Explain your reasoning based on the marker genes
4. List 2-3 alternative possible cell types with their confidence scores

{context_info}

Please respond in the following JSON format:
{{
    "primary_cell_type": "name of most likely cell type",
    "confidence": 0.0-1.0,
    "reasoning": "explanation of why these markers indicate this cell type",
    "alternative_types": [
        {{"cell_type": "alternative 1", "confidence": 0.0-1.0}},
        {{"cell_type": "alternative 2", "confidence": 0.0-1.0}}
    ]
}}

Respond ONLY with valid JSON, no additional text."""
        
        return prompt
    
    def predict_cell_type(
        self,
        marker_genes: List[str],
        cluster_id: str,
        species: str = "human",
        tissue: Optional[str] = None,
        context: Optional[str] = None
    ) -> CellTypePrediction:
        """
        Predict cell type for a cluster based on marker genes.
        
        Parameters
        ----------
        marker_genes : list of str
            Top marker genes for the cluster.
        cluster_id : str
            Cluster identifier.
        species : str
            Species (human, mouse, etc.).
        tissue : str, optional
            Tissue type if known.
        context : str, optional
            Additional context about the experiment.
        
        Returns
        -------
        CellTypePrediction
            Prediction results with confidence and reasoning.
        """
        prompt = self._build_cell_type_prompt(
            marker_genes, cluster_id, species, tissue, context
        )
        
        # API call with retries
        for attempt in range(self.max_retries):
            try:
                self._rate_limit()
                
                response = self.client.chat.completions.create(
                    model=self.model,
                    messages=[
                        {"role": "system", "content": "You are an expert in single-cell biology and cell type identification."},
                        {"role": "user", "content": prompt}
                    ],
                    temperature=self.temperature,
                    max_tokens=500
                )
                
                # Parse response
                content = response.choices[0].message.content.strip()
                
                # Handle potential markdown code blocks
                if content.startswith("```"):
                    content = content.split("```")[1]
                    if content.startswith("json"):
                        content = content[4:]
                    content = content.strip()
                
                result = json.loads(content)
                
                # Build prediction object
                alternative_types = [
                    (alt["cell_type"], alt["confidence"])
                    for alt in result.get("alternative_types", [])
                ]
                
                prediction = CellTypePrediction(
                    cluster_id=cluster_id,
                    predicted_type=result["primary_cell_type"],
                    confidence=float(result["confidence"]),
                    reasoning=result["reasoning"],
                    marker_genes_used=marker_genes,
                    alternative_types=alternative_types
                )
                
                return prediction
                
            except json.JSONDecodeError as e:
                print(f"Warning: Failed to parse JSON response (attempt {attempt + 1}/{self.max_retries})")
                print(f"Response was: {content[:200]}...")
                if attempt < self.max_retries - 1:
                    time.sleep(self.retry_delay)
                else:
                    raise ValueError(f"Failed to get valid JSON response after {self.max_retries} attempts") from e
            
            except Exception as e:
                print(f"Warning: API call failed (attempt {attempt + 1}/{self.max_retries}): {e}")
                if attempt < self.max_retries - 1:
                    time.sleep(self.retry_delay)
                else:
                    raise
    
    def batch_predict_cell_types(
        self,
        clusters: Dict[str, List[str]],
        species: str = "human",
        tissue: Optional[str] = None,
        context: Optional[str] = None,
        progress_callback: Optional[callable] = None
    ) -> Dict[str, CellTypePrediction]:
        """
        Predict cell types for multiple clusters.
        
        Parameters
        ----------
        clusters : dict
            Dictionary mapping cluster_id -> list of marker genes.
        species : str
            Species.
        tissue : str, optional
            Tissue type.
        context : str, optional
            Additional context.
        progress_callback : callable, optional
            Function to call with progress updates: callback(current, total).
        
        Returns
        -------
        dict
            Dictionary mapping cluster_id -> CellTypePrediction.
        """
        predictions = {}
        total = len(clusters)
        
        for idx, (cluster_id, marker_genes) in enumerate(clusters.items(), 1):
            print(f"Processing cluster {cluster_id} ({idx}/{total})...")
            
            prediction = self.predict_cell_type(
                marker_genes=marker_genes,
                cluster_id=cluster_id,
                species=species,
                tissue=tissue,
                context=context
            )
            
            predictions[cluster_id] = prediction
            
            if progress_callback:
                progress_callback(idx, total)
        
        return predictions

# scRN_AI v0.1.0
# Any usage is subject to this software's license.
